// This code is part of the project "Glign: Taming Misaligned Graph Traversals in Concurrent Graph Processing"
// Copyright (c) 2022 Xizhe Yin

// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#define WEIGHTED 1
#include "ligra.h"

using IdxType = long long;

struct SSNP_F {
  intE* WidestPathVal;
  bool* CurrActiveArray;
  bool* NextActiveArray;
  long BatchSize;

  SSNP_F(intE* _WidestPathVal, bool* _CurrActiveArray, bool* _NextActiveArray, long _BatchSize) : 
    WidestPathVal(_WidestPathVal), CurrActiveArray(_CurrActiveArray), NextActiveArray(_NextActiveArray), BatchSize(_BatchSize) {}
  
  inline bool update (uintE s, uintE d, intE edgeLen) { //Update
    bool ret = false;
    IdxType s_begin = s * BatchSize;
    IdxType d_begin = d * BatchSize;

    for (long j = 0; j < BatchSize; j++) {
      if (CurrActiveArray[s_begin + j]) {
        intE sValue = WidestPathVal[s_begin + j];
        intE dValue = WidestPathVal[d_begin + j];
        intE newValue = std::max(sValue, edgeLen);
        if (dValue > newValue) {
          // Visited[d_begin + j] = true;
          WidestPathVal[d_begin + j] = newValue;
          NextActiveArray[d_begin + j] = true;
          ret = true;
        }
      }
    }
    return ret;
  }
  
  inline bool updateAtomic (uintE s, uintE d, intE edgeLen){ //atomic version of Update
    bool ret = false;
    IdxType s_begin = s * BatchSize;
    IdxType d_begin = d * BatchSize;
    for (long j = 0; j < BatchSize; j++) {
      if (CurrActiveArray[s_begin + j]) {
        intE sValue = WidestPathVal[s_begin + j];
        intE newValue = std::max(sValue, edgeLen);
        if (writeMin(&WidestPathVal[d_begin + j], newValue) && CAS(&NextActiveArray[d_begin + j], false, true)) {
          ret = true;
        }
      }
    }
    return ret;
  }
  //cond function checks if vertex has been visited yet
  inline bool cond (uintE d) { return cond_true(d); } 
};

struct SSNP_SKIP_F {
  intE* WidestPathVal;
  bool* CurrActiveArray;
  bool* NextActiveArray;
  long BatchSize;

  SSNP_SKIP_F(intE* _WidestPathVal, long _BatchSize) : 
    WidestPathVal(_WidestPathVal), BatchSize(_BatchSize) {}
  
  inline bool update (uintE s, uintE d, intE edgeLen) { //Update
    bool ret = false;
    IdxType s_begin = s * BatchSize;
    IdxType d_begin = d * BatchSize;

    for (long j = 0; j < BatchSize; j++) {
      intE sValue = WidestPathVal[s_begin + j];
      intE dValue = WidestPathVal[d_begin + j];
      intE newValue = std::max(sValue, edgeLen);
      if (dValue > newValue) {
        // Visited[d_begin + j] = true;
        WidestPathVal[d_begin + j] = newValue;
        ret = true;
      }
    }
    return ret;
  }
  
  inline bool updateAtomic (uintE s, uintE d, intE edgeLen){ //atomic version of Update
    bool ret = false;
    IdxType s_begin = s * BatchSize;
    IdxType d_begin = d * BatchSize;
    for (long j = 0; j < BatchSize; j++) {
      intE sValue = WidestPathVal[s_begin + j];
      intE newValue = std::max(sValue, edgeLen);
      if (writeMin(&WidestPathVal[d_begin + j], newValue)) {
        ret = true;
      }
    }
    return ret;
  }
  //cond function checks if vertex has been visited yet
  inline bool cond (uintE d) { return cond_true(d); } 
};

struct BFSLV_F {
  uintE* Levels;
  bool* CurrActiveArray;
  bool* NextActiveArray;
  long BatchSize;

  BFSLV_F(uintE* _Levels, bool* _CurrActiveArray, bool* _NextActiveArray, long _BatchSize) : 
    Levels(_Levels), CurrActiveArray(_CurrActiveArray), NextActiveArray(_NextActiveArray), BatchSize(_BatchSize) {}
  
  inline bool update (uintE s, uintE d, intE edgeLen=1) { //Update
    bool ret = false;
    IdxType s_begin = s * BatchSize;
    IdxType d_begin = d * BatchSize;

    for (long j = 0; j < BatchSize; j++) {
      if (CurrActiveArray[s_begin + j]) {
        uintE newValue = Levels[s_begin + j] + 1;
        if (Levels[d_begin + j] > newValue) {
          // Visited[d_begin + j] = true;
          Levels[d_begin + j] = newValue;
          NextActiveArray[d_begin + j] = true;
          ret = true;
        }
      }

    }
    return ret;
  }
  
  inline bool updateAtomic (uintE s, uintE d, intE edgeLen=1){ //atomic version of Update
    bool ret = false;
    IdxType s_begin = s * BatchSize;
    IdxType d_begin = d * BatchSize;
    for (long j = 0; j < BatchSize; j++) {
      if (CurrActiveArray[s_begin + j]) {
        uintE newValue = Levels[s_begin + j] + 1;
        if (writeMin(&Levels[d_begin + j], newValue) && CAS(&NextActiveArray[d_begin + j], false, true)) {
          ret = true;
        }
      }
    }
    return ret;
  }
  //cond function checks if vertex has been visited yet
  inline bool cond (uintE d) { return cond_true(d); } 
};

template <class vertex>
uintE* Compute_Eval(graph<vertex>& G, std::vector<long> vecQueries, commandLine P) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size;
  uintE* Hops = pbbs::new_array<uintE>(totalNumVertices);
  bool* CurrActiveArray = pbbs::new_array<bool>(totalNumVertices);
  bool* NextActiveArray = pbbs::new_array<bool>(totalNumVertices);
  bool* frontier = pbbs::new_array<bool>(n);
  parallel_for(size_t i = 0; i < n; i++) {
    frontier[i] = false;
  }
  for(long i = 0; i < batch_size; i++) {
    frontier[vecQueries[i]] = true;
  }
  parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
    Hops[i] = (uintE)MAXLEVEL;
    CurrActiveArray[i] = false;
    NextActiveArray[i] = false;
  }
  for(long i = 0; i < batch_size; i++) {
    Hops[(IdxType)batch_size * (IdxType)vecQueries[i] + (IdxType)i] = 0;
  }
  parallel_for(size_t i = 0; i < batch_size; i++) {
    CurrActiveArray[(IdxType)vecQueries[i] * (IdxType)batch_size + (IdxType)i] = true;
  }

  vertexSubset Frontier(n, frontier);
  while(!Frontier.isEmpty()){
    // mode: no_dense, remove_duplicates (for batch size > 1)
    vertexSubset output = edgeMap(G, Frontier, BFSLV_F(Hops, CurrActiveArray, NextActiveArray, batch_size), -1, no_dense|remove_duplicates);

    Frontier.del();
    Frontier = output;

    std::swap(CurrActiveArray, NextActiveArray);
    parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
      NextActiveArray[i] = false;
    }
  }

  Frontier.del();
  pbbs::delete_array(CurrActiveArray, totalNumVertices);
  pbbs::delete_array(NextActiveArray, totalNumVertices);
  return Hops;
}

template <class vertex>
pair<size_t, size_t> Compute_Base(graph<vertex>& G, std::vector<long> vecQueries, commandLine P, bool should_profile) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size;
  intE* WidestPathVal = pbbs::new_array<intE>(totalNumVertices);
  bool* CurrActiveArray = pbbs::new_array<bool>(totalNumVertices);
  bool* NextActiveArray = pbbs::new_array<bool>(totalNumVertices);
  bool* frontier = pbbs::new_array<bool>(n);
  parallel_for(size_t i = 0; i < n; i++) {
    frontier[i] = false;
  }
  for(long i = 0; i < batch_size; i++) {
    frontier[vecQueries[i]] = true;
  }
  parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
    WidestPathVal[i] = (intE)MAXWIDTH;
    CurrActiveArray[i] = false;
    NextActiveArray[i] = false;
  }
  for(long i = 0; i < batch_size; i++) {
    WidestPathVal[(IdxType)batch_size * (IdxType)vecQueries[i] + (IdxType)i] = (intE)0;
  }
  parallel_for(size_t i = 0; i < batch_size; i++) {
    CurrActiveArray[(IdxType)vecQueries[i] * (IdxType)batch_size + (IdxType)i] = true;
  }
  vertexSubset Frontier(n, frontier);
  // for profiling
  long iteration = 0;
  size_t totalActivated = 0;
  while(!Frontier.isEmpty()){
    iteration++;
    totalActivated += Frontier.size();
    // mode: no_dense, remove_duplicates (for batch size > 1)
    vertexSubset output = edgeMap(G, Frontier, SSNP_F(WidestPathVal, CurrActiveArray, NextActiveArray, batch_size), -1, no_dense|remove_duplicates);

    Frontier.del();
    Frontier = output;

    std::swap(CurrActiveArray, NextActiveArray);
    parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
      NextActiveArray[i] = false;
    }
  }

#ifdef OUTPUT 
  for (int i = 0; i < batch_size; i++) {
    long start = vecQueries[i];
    char outFileName[300];
    sprintf(outFileName, "SSNP_base_output_src%ld.%ld.%ld.out", start, edge_count, batch_size);
    FILE *fp;
    fp = fopen(outFileName, "w");
    for (long j = 0; j < n; j++)
      fprintf(fp, "%ld %d\n", j, WidestPathVal[j * batch_size + i]);
    fclose(fp);
  }
#endif

  Frontier.del();
  pbbs::delete_array(WidestPathVal, totalNumVertices);
  pbbs::delete_array(CurrActiveArray, totalNumVertices);
  pbbs::delete_array(NextActiveArray, totalNumVertices);
  return make_pair(totalActivated, 0);
}


template <class vertex>
pair<size_t, size_t> Compute_Base_Skipping(graph<vertex>& G, std::vector<long> vecQueries, commandLine P, bool should_profile) {
    size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size;
  intE* WidestPathVal = pbbs::new_array<intE>(totalNumVertices);
  bool* frontier = pbbs::new_array<bool>(n);
  parallel_for(size_t i = 0; i < n; i++) {
    frontier[i] = false;
  }
  for(long i = 0; i < batch_size; i++) {
    frontier[vecQueries[i]] = true;
  }
  parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
    WidestPathVal[i] = (intE)MAXWIDTH;
  }
  for(long i = 0; i < batch_size; i++) {
    WidestPathVal[(IdxType)batch_size * (IdxType)vecQueries[i] + (IdxType)i] = (intE)0;
  }

  vertexSubset Frontier(n, frontier);
  // for profiling
  long iteration = 0;
  size_t totalActivated = 0;

  while(!Frontier.isEmpty()){
    iteration++;
    totalActivated += Frontier.size();
    // mode: no_dense, remove_duplicates (for batch size > 1)
    vertexSubset output = edgeMap(G, Frontier, SSNP_SKIP_F(WidestPathVal, batch_size), -1, no_dense|remove_duplicates);
    Frontier.del();
    Frontier = output;
  }

#ifdef OUTPUT 
  for (int i = 0; i < batch_size; i++) {
    long start = vecQueries[i];
    char outFileName[300];
    sprintf(outFileName, "SSNP_base_skipping_output_src%ld.%ld.%ld.out", start, edge_count, batch_size);
    FILE *fp;
    fp = fopen(outFileName, "w");
    for (long j = 0; j < n; j++)
      fprintf(fp, "%ld %d\n", j, WidestPathVal[j * batch_size + i]);
    fclose(fp);
  }
#endif

  Frontier.del();
  pbbs::delete_array(WidestPathVal, totalNumVertices);
  return make_pair(totalActivated, 0);
}
template <class vertex>
pair<size_t, size_t> Compute_Delay_Skipping(graph<vertex>& G, std::vector<long> vecQueries, commandLine P, std::vector<int> defer_vec, bool should_profile) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size;
  intE* WidestPathVal = pbbs::new_array<intE>(totalNumVertices);
  bool* frontier = pbbs::new_array<bool>(n);
  parallel_for(size_t i = 0; i < n; i++) {
    frontier[i] = false;
  }
  // for delaying initialization
  for(long i = 0; i < batch_size; i++) {
    if (defer_vec[i] == 0) {
      frontier[vecQueries[i]] = true;
    }
  }

  parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
    WidestPathVal[i] = (intE)MAXWIDTH;
  }
  for(long i = 0; i < batch_size; i++) {
    WidestPathVal[(IdxType)batch_size * (IdxType)vecQueries[i] + (IdxType)i] = (intE)0;
  }

  vertexSubset Frontier(n, frontier);

  // for profiling
  long iteration = 0;
  size_t totalActivated = 0;

  while(!Frontier.isEmpty()){
    iteration++;
    totalActivated += Frontier.size();

    // mode: no_dense, remove_duplicates (for batch size > 1)
    vertexSubset output = edgeMap(G, Frontier, SSNP_SKIP_F(WidestPathVal, batch_size), -1, no_dense|remove_duplicates);

    Frontier.del();
    Frontier = output;

    Frontier.toDense();
    bool* new_d = Frontier.d;
    Frontier.d = nullptr;
    for(long i = 0; i < batch_size; i++) {
      if (defer_vec[i] == iteration) {
        new_d[vecQueries[i]] = true;
      }
    }
    vertexSubset Frontier_new(n, new_d);
    Frontier.del();
    Frontier = Frontier_new;
  }


#ifdef OUTPUT 
  for (int i = 0; i < batch_size; i++) {
    long start = vecQueries[i];
    char outFileName[300];
    sprintf(outFileName, "SSNP_delay_skipping_output_src%ld.%ld.%ld.out", start, edge_count, batch_size);
    FILE *fp;
    fp = fopen(outFileName, "w");
    for (long j = 0; j < n; j++)
      fprintf(fp, "%ld %d\n", j, WidestPathVal[j * batch_size + i]);
    fclose(fp);
  }
#endif

  Frontier.del();
  pbbs::delete_array(WidestPathVal, totalNumVertices);
  return make_pair(totalActivated, 0);
}
template <class vertex>
pair<size_t, size_t> Compute_Base_Dynamic(graph<vertex>& G, std::vector<long> vecQueries, queue<long>& queryQueue, commandLine P, bool should_profile) {
  return make_pair(0,0);
}

template <class vertex>
pair<size_t, size_t> Compute_Delay(graph<vertex>& G, std::vector<long> vecQueries, commandLine P, std::vector<int> defer_vec, bool should_profile) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size;
  intE* WidestPathVal = pbbs::new_array<intE>(totalNumVertices);
  bool* CurrActiveArray = pbbs::new_array<bool>(totalNumVertices);
  bool* NextActiveArray = pbbs::new_array<bool>(totalNumVertices);
  bool* frontier = pbbs::new_array<bool>(n);
  parallel_for(size_t i = 0; i < n; i++) {
    frontier[i] = false;
  }
  // for delaying initialization
  for(long i = 0; i < batch_size; i++) {
    if (defer_vec[i] == 0) {
      frontier[vecQueries[i]] = true;
    }
  }

  parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
    WidestPathVal[i] = (intE)MAXWIDTH;
    CurrActiveArray[i] = false;
    NextActiveArray[i] = false;
  }
  for(long i = 0; i < batch_size; i++) {
    WidestPathVal[(IdxType)batch_size * (IdxType)vecQueries[i] + (IdxType)i] = (intE)0;
  }
  parallel_for(size_t i = 0; i < batch_size; i++) {
    if (defer_vec[i] == 0) {
      CurrActiveArray[(IdxType)vecQueries[i] * (IdxType)batch_size + (IdxType)i] = true;
    }
  }

  vertexSubset Frontier(n, frontier);

  // for profiling
  long iteration = 0;
  size_t totalActivated = 0;

  while(!Frontier.isEmpty()){
    iteration++;
    totalActivated += Frontier.size();

    // mode: no_dense, remove_duplicates (for batch size > 1)
    vertexSubset output = edgeMap(G, Frontier, SSNP_F(WidestPathVal, CurrActiveArray, NextActiveArray, batch_size), -1, no_dense|remove_duplicates);

    Frontier.del();
    Frontier = output;

    Frontier.toDense();
    bool* new_d = Frontier.d;
    Frontier.d = nullptr;
    for(long i = 0; i < batch_size; i++) {
      if (defer_vec[i] == iteration) {
        new_d[vecQueries[i]] = true;
        NextActiveArray[(IdxType)vecQueries[i] * (IdxType)batch_size + (IdxType)i] = true;
      }
    }
    vertexSubset Frontier_new(n, new_d);
    Frontier.del();
    Frontier = Frontier_new;

    std::swap(CurrActiveArray, NextActiveArray);
    parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
      NextActiveArray[i] = false;
    }
  }


#ifdef OUTPUT 
  for (int i = 0; i < batch_size; i++) {
    long start = vecQueries[i];
    char outFileName[300];
    sprintf(outFileName, "SSNP_delay_output_src%ld.%ld.%ld.out", start, edge_count, batch_size);
    FILE *fp;
    fp = fopen(outFileName, "w");
    for (long j = 0; j < n; j++)
      fprintf(fp, "%ld %d\n", j, WidestPathVal[j * batch_size + i]);
    fclose(fp);
  }
#endif

  Frontier.del();
  pbbs::delete_array(CurrActiveArray, totalNumVertices);
  pbbs::delete_array(NextActiveArray, totalNumVertices);
  pbbs::delete_array(WidestPathVal, totalNumVertices);
  return make_pair(totalActivated, 0);
}