// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of 
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
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
// typedef uint32_t benchmarkType;
// const benchmarkType _empty = 0;
// const benchmarkType benchmark_sssp = 1;
// const benchmarkType benchmark_sswp = 2;
// const benchmarkType benchmark_viterbi = 4;
// const benchmarkType benchmark_bfs = 8;

struct DJ_F {
  intE* ShortestPathLen;
  bool* CurrActiveArray;
  bool* NextActiveArray;
  long BatchSize;

  DJ_F(intE* _ShortestPathLen, bool* _CurrActiveArray, bool* _NextActiveArray, long _BatchSize) : 
    ShortestPathLen(_ShortestPathLen), CurrActiveArray(_CurrActiveArray), NextActiveArray(_NextActiveArray), BatchSize(_BatchSize) {}
  
  inline bool update (uintE s, uintE d, intE edgeLen) { //Update
    bool ret = false;
    IdxType s_begin = s * BatchSize;
    IdxType d_begin = d * BatchSize;

    for (long j = 0; j < BatchSize; j++) {
      if (CurrActiveArray[s_begin + j]) {
        intE newValue = ShortestPathLen[s_begin + j] + edgeLen;
        if (ShortestPathLen[d_begin + j] > newValue) {
          // Visited[d_begin + j] = true;
          ShortestPathLen[d_begin + j] = newValue;
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
        intE newValue = ShortestPathLen[s_begin + j] + edgeLen;
        if (writeMin(&ShortestPathLen[d_begin + j], newValue) && CAS(&NextActiveArray[d_begin + j], false, true)) {
          ret = true;
        }
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


struct Heterogeneous_F {
    intE* QueryValues;
    benchmarkType* CurrActiveArray;
    benchmarkType* NextActiveArray;
    long BatchSize;

    Heterogeneous_F(intE* _QueryValues, benchmarkType* _CurrActiveArray, benchmarkType* _NextActiveArray, long _BatchSize) : 
        QueryValues(_QueryValues), CurrActiveArray(_CurrActiveArray), NextActiveArray(_NextActiveArray), BatchSize(_BatchSize) {}

    inline bool update(uintE s, uintE d, intE edgeLen=1) {
        bool ret = false;
        IdxType s_begin = s * BatchSize;
        IdxType d_begin = d * BatchSize;
        for (long j = 0; j < BatchSize; j++) {
          if (CurrActiveArray[s_begin + j] == benchmark_sssp) {
            intE sValue = QueryValues[s_begin + j];
            intE dValue = QueryValues[d_begin + j];
            intE newDist = sValue + edgeLen;
            if (dValue > newDist) {
                QueryValues[d_begin + j] = newDist;
                NextActiveArray[d_begin + j] = benchmark_sssp;
                ret = true;
            }
            continue;
          } else if (CurrActiveArray[s_begin + j] == benchmark_bfs) {
            intE sValue = QueryValues[s_begin + j];
            intE dValue = QueryValues[d_begin + j];
            intE newDist = sValue + 1;
            if (dValue > newDist) {
                QueryValues[d_begin + j] = newDist;
                NextActiveArray[d_begin + j] = benchmark_bfs;
                ret = true;
            }
            continue;
          } else if (CurrActiveArray[s_begin + j] == benchmark_sswp) {
            intE sValue = QueryValues[s_begin + j];
            intE dValue = QueryValues[d_begin + j];
            intE newWidth = std::min(sValue, edgeLen);
            if (dValue < newWidth) {
                QueryValues[d_begin + j] = newWidth;
                NextActiveArray[d_begin + j] = benchmark_sswp;
                ret = true;
            }
            continue;
          } else if (CurrActiveArray[s_begin + j] == benchmark_ssnp) {
            intE sValue = QueryValues[s_begin + j];
            intE dValue = QueryValues[d_begin + j];
            intE newWidth = std::max(sValue, edgeLen);
            if (dValue > newWidth) {
                QueryValues[d_begin + j] = newWidth;
                NextActiveArray[d_begin + j] = benchmark_ssnp;
                ret = true;
            }
            continue;
          }
        }
        return ret;
    }

    inline bool updateAtomic(uintE s, uintE d, intE edgeLen=1) {
      bool ret = false;
      IdxType s_begin = s * BatchSize;
      IdxType d_begin = d * BatchSize;
      for (long j = 0; j < BatchSize; j++) {
        auto cur_active = CurrActiveArray[s_begin + j];
        if (cur_active == benchmark_sssp) {
          intE sValue = QueryValues[s_begin + j];
          intE newDist = sValue + edgeLen;
          if (writeMin(&QueryValues[d_begin + j], newDist)
                  && CAS(&NextActiveArray[d_begin + j], _empty, benchmark_sssp)) {
              ret = true;
          }
          continue;
        } else if (cur_active == benchmark_bfs) {
          intE sValue = QueryValues[s_begin + j];
          intE newDist = sValue + 1;
          if (writeMin(&QueryValues[d_begin + j], newDist)
                  && CAS(&NextActiveArray[d_begin + j], _empty, benchmark_bfs)) {
              ret = true;
          }
          continue;
        } else if (cur_active == benchmark_sswp) {
          intE sValue = QueryValues[s_begin + j];
          intE newWidth = std::min(sValue, edgeLen);
          if (writeMax(&QueryValues[d_begin + j], newWidth)
                  && CAS(&NextActiveArray[d_begin + j], _empty, benchmark_sswp)) {
              ret = true;
          }
          continue;
        } else if (cur_active == benchmark_ssnp) {
          intE sValue = QueryValues[s_begin + j];
          intE newWidth = std::max(sValue, edgeLen);
          if (writeMin(&QueryValues[d_begin + j], newWidth)
                  && CAS(&NextActiveArray[d_begin + j], _empty, benchmark_ssnp)) {
              ret = true;
          }
          continue;
        }
      }
      return ret;
    }

    inline bool cond_true (int d) { return 1; }

    inline bool cond(uintE d) {
      return cond_true(d);
    }
};

struct Heterogeneous_SKIP_F {
    intE* QueryValues;
    long BatchSize;
    benchmarkType* typeArray;

    Heterogeneous_SKIP_F(intE* _QueryValues, long _BatchSize, benchmarkType* _typeArray) : 
      QueryValues(_QueryValues), BatchSize(_BatchSize), typeArray(_typeArray) {}

    inline bool update(uintE s, uintE d, intE edgeLen=1) {
      bool ret = false;
      IdxType s_begin = s * BatchSize;
      IdxType d_begin = d * BatchSize;
      for (long j = 0; j < BatchSize; j++) {
        if (typeArray[j] == benchmark_sssp) 
        {
          intE sValue = QueryValues[s_begin + j];
          intE dValue = QueryValues[d_begin + j];
          intE newDist = sValue + edgeLen;
          if (dValue > newDist) {
              QueryValues[d_begin + j] = newDist;
              ret = true;
          }
          continue;
        } 
        else if (typeArray[j] == benchmark_bfs) 
        {
          intE sValue = QueryValues[s_begin + j];
          intE dValue = QueryValues[d_begin + j];
          intE newDist = sValue + 1;

          if (dValue > newDist) {
              QueryValues[d_begin + j] = newDist;
              ret = true;
          }
          continue;
        } 
        else if (typeArray[j] == benchmark_sswp) 
        {
          intE sValue = QueryValues[s_begin + j];
          intE dValue = QueryValues[d_begin + j];
          intE newWidth = std::min(sValue, edgeLen);
          if (dValue < newWidth) {
              QueryValues[d_begin + j] = newWidth;
              ret = true;                    
          }
          continue;
        }  
        else if (typeArray[j] == benchmark_ssnp) 
        {
          intE sValue = QueryValues[s_begin + j];
          intE dValue = QueryValues[d_begin + j];
          intE newWidth = std::max(sValue, edgeLen);
          if (dValue > newWidth) {
              QueryValues[d_begin + j] = newWidth;
              ret = true;                    
          }
          continue;
        } 
      }
      return ret;
    }

    inline bool updateAtomic(uintE s, uintE d, intE edgeLen=1) {
      bool ret = false;
      IdxType s_begin = s * BatchSize;
      IdxType d_begin = d * BatchSize;
      for (long j = 0; j < BatchSize; j++) {
        if (typeArray[j] == benchmark_sssp) 
        {
          intE sValue = QueryValues[s_begin + j];
          intE newDist = sValue + edgeLen;
          if (writeMin(&QueryValues[d_begin + j], newDist)) {
              ret = true;
          }
          continue;
        } 
        else if (typeArray[j] == benchmark_bfs) 
        {
          intE sValue = QueryValues[s_begin + j];
          intE newDist = sValue + 1;
          if (writeMin(&QueryValues[d_begin + j], newDist)) {
              ret = true;
          }
          continue;
        } 
        else if (typeArray[j] == benchmark_sswp) 
        {
          intE sValue = QueryValues[s_begin + j];
          intE newWidth = std::min(sValue, edgeLen);
          if (writeMax(&QueryValues[d_begin + j], newWidth)) {
              ret = true;
          }
          continue;
        } 
        else if (typeArray[j] == benchmark_ssnp) 
        {
          intE sValue = QueryValues[s_begin + j];
          intE newWidth = std::max(sValue, edgeLen);
          if (writeMin(&QueryValues[d_begin + j], newWidth)) {
              ret = true;
          }
          continue;
        } 
      }
      return ret;
    }

    inline bool cond_true (int d) { return 1; }

    inline bool cond(uintE d) {
      return cond_true(d);
    }
};

template <class vertex>
uintE* Compute_Eval(graph<vertex>& G, std::vector<long> vecQueries, commandLine P) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size;
  uintE* Levels = pbbs::new_array<uintE>(totalNumVertices);
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
    Levels[i] = (uintE)MAXLEVEL;
    CurrActiveArray[i] = false;
    NextActiveArray[i] = false;
  }
  for(long i = 0; i < batch_size; i++) {
    Levels[(IdxType)batch_size * (IdxType)vecQueries[i] + (IdxType)i] = 0;
  }
  parallel_for(size_t i = 0; i < batch_size; i++) {
    CurrActiveArray[(IdxType)vecQueries[i] * (IdxType)batch_size + (IdxType)i] = true;
  }

  vertexSubset Frontier(n, frontier);
  while(!Frontier.isEmpty()){
    // mode: no_dense, remove_duplicates (for batch size > 1)
    vertexSubset output = edgeMap(G, Frontier, BFSLV_F(Levels, CurrActiveArray, NextActiveArray, batch_size), -1, no_dense|remove_duplicates);

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
  return Levels;
}


template <class vertex>
pair<size_t, size_t> Compute_Heter_Base(graph<vertex>& G, std::vector<pair<long, benchmarkType>> vecQueries, commandLine P, bool should_profile) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size;
  intE* QueryValues = pbbs::new_array<intE>(totalNumVertices);
  benchmarkType* CurrActiveArray = pbbs::new_array<benchmarkType>(totalNumVertices);
  benchmarkType* NextActiveArray = pbbs::new_array<benchmarkType>(totalNumVertices);
  bool* frontier = pbbs::new_array<bool>(n);
  parallel_for(size_t i = 0; i < n; i++) {
    frontier[i] = false;
  }
  for(long i = 0; i < batch_size; i++) {
    frontier[vecQueries[i].first] = true;
  }
  parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
    CurrActiveArray[i] = _empty;
    NextActiveArray[i] = _empty;
  }
  parallel_for(IdxType i = 0; i < batch_size; i++) {
    if (vecQueries[i].second & benchmark_sssp) {
      CurrActiveArray[(IdxType)vecQueries[i].first * (IdxType)batch_size + (IdxType)i] = benchmark_sssp;
    } else if (vecQueries[i].second & benchmark_bfs) {
      CurrActiveArray[(IdxType)vecQueries[i].first * (IdxType)batch_size + (IdxType)i] = benchmark_bfs;
    } else if (vecQueries[i].second & benchmark_sswp) {
      CurrActiveArray[(IdxType)vecQueries[i].first * (IdxType)batch_size + (IdxType)i] = benchmark_sswp;
    } else if (vecQueries[i].second & benchmark_ssnp) {
      CurrActiveArray[(IdxType)vecQueries[i].first * (IdxType)batch_size + (IdxType)i] = benchmark_ssnp;
    }
  }
  for(long i = 0; i < batch_size; i++) {
    if (vecQueries[i].second & benchmark_sssp) {
      parallel_for(size_t j = 0; j < n; j++) {
        QueryValues[j * batch_size + i] = (intE)MAXPATH;
      }
    } else if (vecQueries[i].second & benchmark_bfs) {
      parallel_for(size_t j = 0; j < n; j++) {
        QueryValues[j * batch_size + i] = (uintE)MAXLEVEL;
      }
    } else if (vecQueries[i].second & benchmark_sswp) {
      parallel_for(size_t j = 0; j < n; j++) {
        QueryValues[j * batch_size + i] = (intE)0;
      }
    } else if (vecQueries[i].second & benchmark_ssnp) {
      parallel_for(size_t j = 0; j < n; j++) {
        QueryValues[j * batch_size + i] = (intE)MAXWIDTH;
      }
    } 
  }

  for(size_t i = 0; i < batch_size; i++) {
    if (vecQueries[i].second & benchmark_sssp) {
      QueryValues[(IdxType)batch_size * (IdxType)vecQueries[i].first + (IdxType)i] = 0;
    } else if (vecQueries[i].second & benchmark_bfs) {
      QueryValues[(IdxType)batch_size * (IdxType)vecQueries[i].first + (IdxType)i] = 0;
    } else if (vecQueries[i].second & benchmark_sswp) {
      QueryValues[(IdxType)batch_size * (IdxType)vecQueries[i].first + (IdxType)i] = (intE)(MAXWIDTH);
    } else if (vecQueries[i].second & benchmark_ssnp) {
      QueryValues[(IdxType)batch_size * (IdxType)vecQueries[i].first + (IdxType)i] = (intE)0;
    }
  }

  vertexSubset Frontier(n, frontier);

  // for profiling
  long iteration = 0;
  size_t totalActivated = 0;
  size_t totalNoOverlap = 0;

  while(!Frontier.isEmpty()){
    iteration++;
    totalActivated += Frontier.size();

    // mode: no_dense, remove_duplicates (for batch size > 1)
    vertexSubset output = edgeMap(G, Frontier, Heterogeneous_F(QueryValues, CurrActiveArray, NextActiveArray, batch_size), -1, no_dense|remove_duplicates);

    Frontier.del();
    Frontier = output;

    Frontier.toDense();
    // bool* new_d = Frontier.d;
    // Frontier.d = nullptr;
    // vertexSubset Frontier_new(n, new_d);
    // Frontier.del();
    // Frontier = Frontier_new;

    std::swap(CurrActiveArray, NextActiveArray);
    parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
      NextActiveArray[i] = false;
    }
  }

  // profiling
  
#ifdef OUTPUT 
  for (int i = 0; i < batch_size; i++) {
    long start = vecQueries[i].first;
    char outFileName[300];
    sprintf(outFileName, "Heter_base_output_src%ld.%ld.%ld.out", start, edge_count, batch_size);
    FILE *fp;
    fp = fopen(outFileName, "w");
    for (long j = 0; j < n; j++)
      fprintf(fp, "%ld %d\n", j, QueryValues[j * batch_size + i]);
    fclose(fp);
  }
#endif

  Frontier.del();
  pbbs::delete_array(QueryValues, totalNumVertices);
  pbbs::delete_array(CurrActiveArray, totalNumVertices);
  pbbs::delete_array(NextActiveArray, totalNumVertices);

  return make_pair(totalActivated, totalNoOverlap);
}

template <class vertex>
pair<size_t, size_t> Compute_Heter_Skip(graph<vertex>& G, std::vector<pair<long, benchmarkType>> vecQueries, commandLine P, bool should_profile) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size;
  intE* QueryValues = pbbs::new_array<intE>(totalNumVertices);
  // benchmarkType* CurrActiveArray = pbbs::new_array<benchmarkType>(totalNumVertices);
  // benchmarkType* NextActiveArray = pbbs::new_array<benchmarkType>(totalNumVertices);
  bool* frontier = pbbs::new_array<bool>(n);
  benchmarkType* typeArray = pbbs::new_array<benchmarkType>(batch_size);
  parallel_for(size_t i = 0; i < n; i++) {
    frontier[i] = false;
  }
  for(long i = 0; i < batch_size; i++) {
    frontier[vecQueries[i].first] = true;
    typeArray[i] = vecQueries[i].second;
  }

  for(long i = 0; i < batch_size; i++) {
    if (vecQueries[i].second & benchmark_sssp) {
      parallel_for(size_t j = 0; j < n; j++) {
        QueryValues[j * batch_size + i] = (intE)MAXPATH;
      }
    } else if (vecQueries[i].second & benchmark_bfs) {
      parallel_for(size_t j = 0; j < n; j++) {
        QueryValues[j * batch_size + i] = (uintE)MAXLEVEL;
      }
    } else if (vecQueries[i].second & benchmark_sswp) {
      parallel_for(size_t j = 0; j < n; j++) {
        QueryValues[j * batch_size + i] = (intE)0;
      }
    } else if (vecQueries[i].second & benchmark_ssnp) {
      parallel_for(size_t j = 0; j < n; j++) {
        QueryValues[j * batch_size + i] = (intE)MAXWIDTH;
      }
    } 
  }

  for(size_t i = 0; i < batch_size; i++) {
    if (vecQueries[i].second & benchmark_sssp) {
      QueryValues[(IdxType)batch_size * (IdxType)vecQueries[i].first + (IdxType)i] = 0;
    } else if (vecQueries[i].second & benchmark_bfs) {
      QueryValues[(IdxType)batch_size * (IdxType)vecQueries[i].first + (IdxType)i] = 0;
    } else if (vecQueries[i].second & benchmark_sswp) {
      QueryValues[(IdxType)batch_size * (IdxType)vecQueries[i].first + (IdxType)i] = (intE)(MAXWIDTH);
    } else if (vecQueries[i].second & benchmark_ssnp) {
      QueryValues[(IdxType)batch_size * (IdxType)vecQueries[i].first + (IdxType)i] = (intE)0;
    }
  }

  vertexSubset Frontier(n, frontier);

  // for profiling
  long iteration = 0;
  size_t totalActivated = 0;
  size_t totalNoOverlap = 0;

  while(!Frontier.isEmpty()){
    iteration++;
    totalActivated += Frontier.size();
    // mode: no_dense, remove_duplicates (for batch size > 1)
    vertexSubset output = edgeMap(G, Frontier, Heterogeneous_SKIP_F(QueryValues, batch_size, typeArray), -1, no_dense|remove_duplicates);

    Frontier.del();
    Frontier = output;

    Frontier.toDense();
  }  
#ifdef OUTPUT 
  for (int i = 0; i < batch_size; i++) {
    long start = vecQueries[i].first;
    char outFileName[300];
    sprintf(outFileName, "Heter_base_skip_output_src%ld.%ld.%ld.out", start, edge_count, batch_size);
    FILE *fp;
    fp = fopen(outFileName, "w");
    for (long j = 0; j < n; j++)
      fprintf(fp, "%ld %d\n", j, QueryValues[j * batch_size + i]);
    fclose(fp);
  }
#endif

  Frontier.del();
  pbbs::delete_array(QueryValues, totalNumVertices);
  pbbs::delete_array(typeArray, batch_size);

  return make_pair(totalActivated, totalNoOverlap);
}

template <class vertex>
pair<size_t, size_t> Compute_Heter_Delay(graph<vertex>& G, std::vector<pair<long, benchmarkType>> vecQueries, commandLine P, std::vector<int> defer_vec, bool should_profile) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size;
  intE* QueryValues = pbbs::new_array<intE>(totalNumVertices);
  bool* frontier = pbbs::new_array<bool>(n);
  benchmarkType* typeArray = pbbs::new_array<benchmarkType>(batch_size);
  parallel_for(size_t i = 0; i < n; i++) {
    frontier[i] = false;
  }
  // for delaying initialization
  for(long i = 0; i < batch_size; i++) {
    typeArray[i] = vecQueries[i].second;
    if (defer_vec[i] == 0) {
      frontier[vecQueries[i].first] = true;
    }
  }

  for(long i = 0; i < batch_size; i++) {
    if (vecQueries[i].second & benchmark_sssp) {
      parallel_for(size_t j = 0; j < n; j++) {
        QueryValues[j * batch_size + i] = (intE)MAXPATH;
      }
    } else if (vecQueries[i].second & benchmark_bfs) {
      parallel_for(size_t j = 0; j < n; j++) {
        QueryValues[j * batch_size + i] = (uintE)MAXLEVEL;
      }
    } else if (vecQueries[i].second & benchmark_sswp) {
      parallel_for(size_t j = 0; j < n; j++) {
        QueryValues[j * batch_size + i] = (intE)0;
      }
    } else if (vecQueries[i].second & benchmark_ssnp) {
      parallel_for(size_t j = 0; j < n; j++) {
        QueryValues[j * batch_size + i] = (intE)MAXWIDTH;
      }
    } 
  }

  for(size_t i = 0; i < batch_size; i++) {
    if (vecQueries[i].second & benchmark_sssp) {
      QueryValues[(IdxType)batch_size * (IdxType)vecQueries[i].first + (IdxType)i] = 0;
    } else if (vecQueries[i].second & benchmark_bfs) {
      QueryValues[(IdxType)batch_size * (IdxType)vecQueries[i].first + (IdxType)i] = 0;
    } else if (vecQueries[i].second & benchmark_sswp) {
      QueryValues[(IdxType)batch_size * (IdxType)vecQueries[i].first + (IdxType)i] = (intE)(MAXWIDTH);
    } else if (vecQueries[i].second & benchmark_ssnp) {
      QueryValues[(IdxType)batch_size * (IdxType)vecQueries[i].first + (IdxType)i] = (intE)0;
    }
  }

  vertexSubset Frontier(n, frontier);

  // for profiling
  long iteration = 0;
  size_t totalActivated = 0;
  size_t totalNoOverlap = 0;

  while(!Frontier.isEmpty()){
    iteration++;
    totalActivated += Frontier.size();

    // mode: no_dense, remove_duplicates (for batch size > 1)
    vertexSubset output = edgeMap(G, Frontier, Heterogeneous_SKIP_F(QueryValues, batch_size, typeArray), -1, no_dense|remove_duplicates);

    Frontier.del();
    Frontier = output;

    Frontier.toDense();
    bool* new_d = Frontier.d;
    Frontier.d = nullptr;
    for(long i = 0; i < batch_size; i++) {
      if (defer_vec[i] == iteration) {
        new_d[vecQueries[i].first] = true;
      }
    }
    vertexSubset Frontier_new(n, new_d);
    Frontier.del();
    Frontier = Frontier_new;
  }

#ifdef OUTPUT 
  for (int i = 0; i < batch_size; i++) {
    long start = vecQueries[i].first;
    char outFileName[300];
    sprintf(outFileName, "Heter_delay_skip_output_src%ld.%ld.%ld.out", start, edge_count, batch_size);
    FILE *fp;
    fp = fopen(outFileName, "w");
    for (long j = 0; j < n; j++)
      fprintf(fp, "%ld %d\n", j, QueryValues[j * batch_size + i]);
    fclose(fp);
  }
#endif

  Frontier.del();
  pbbs::delete_array(QueryValues, totalNumVertices);
  return make_pair(totalActivated, totalNoOverlap);
}

template<class vertex>
void Compute_Delay(graph<vertex>&, vector<long>, commandLine, vector<int>, bool should_profile){

  return;
}