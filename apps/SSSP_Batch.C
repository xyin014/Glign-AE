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

struct DJ_SINGLE_F {
  intE* ShortestPathLen;

  DJ_SINGLE_F(intE* _ShortestPathLen) : 
    ShortestPathLen(_ShortestPathLen) {}
  
  inline bool update (uintE s, uintE d, intE edgeLen) { //Update
    bool ret = false;
    intE newValue = ShortestPathLen[s] + edgeLen;
    if (ShortestPathLen[d] > newValue) {
      ShortestPathLen[d] = newValue;
      ret = true;
    }
    return ret;
  }
  
  inline bool updateAtomic (uintE s, uintE d, intE edgeLen){ //atomic version of Update
    bool ret = false;
    IdxType s_begin = s;
    IdxType d_begin = d;
    intE newValue = ShortestPathLen[s] + edgeLen;
    if (writeMin(&ShortestPathLen[d], newValue)) {
      ret = true;
    }
    return ret;
  }
  //cond function checks if vertex has been visited yet
  inline bool cond (uintE d) { return cond_true(d); } 
};

struct DJ_SKIP_F {
  intE* ShortestPathLen;
  long BatchSize;

  DJ_SKIP_F(intE* _ShortestPathLen, long _BatchSize) : 
    ShortestPathLen(_ShortestPathLen), BatchSize(_BatchSize) {}
  
  inline bool update (uintE s, uintE d, intE edgeLen) { //Update
    bool ret = false;
    IdxType s_begin = s * BatchSize;
    IdxType d_begin = d * BatchSize;

    for (long j = 0; j < BatchSize; j++) {
      // if (CurrActiveArray[s_begin + j]) {
        intE newValue = ShortestPathLen[s_begin + j] + edgeLen;
        if (ShortestPathLen[d_begin + j] > newValue) {
          // Visited[d_begin + j] = true;
          ShortestPathLen[d_begin + j] = newValue;
          // NextActiveArray[d_begin + j] = true;
          ret = true;
        }
      // }
    }
    return ret;
  }
  
  inline bool updateAtomic (uintE s, uintE d, intE edgeLen){ //atomic version of Update
    bool ret = false;
    IdxType s_begin = s * BatchSize;
    IdxType d_begin = d * BatchSize;
    for (long j = 0; j < BatchSize; j++) {
      // if (CurrActiveArray[s_begin + j]) {
        intE newValue = ShortestPathLen[s_begin + j] + edgeLen;
        if (writeMin(&ShortestPathLen[d_begin + j], newValue) ) {
          ret = true;
        }
      // }
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
uintE* Compute_Eval_Prop(graph<vertex>& G, std::vector<long> vecQueries, commandLine P) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size;
  intE* Levels = pbbs::new_array<intE>(totalNumVertices);
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
    Levels[i] = (intE)MAXLEVEL;
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
    vertexSubset output = edgeMap(G, Frontier, DJ_F(Levels, CurrActiveArray, NextActiveArray, batch_size), -1, no_dense|remove_duplicates);

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
  uintE* ret = pbbs::new_array<uintE>(totalNumVertices);
  parallel_for(size_t i = 0; i < totalNumVertices; i++) {
    ret[i] = (uintE)Levels[i];
  }
  pbbs::delete_array(Levels, totalNumVertices);
  return ret;
}

template <class vertex>
pair<size_t, size_t> Compute_Base(graph<vertex>& G, std::vector<long> vecQueries, commandLine P, bool should_profile) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size;
  intE* ShortestPathLen = pbbs::new_array<intE>(totalNumVertices);
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
    ShortestPathLen[i] = (intE)MAXPATH;
    CurrActiveArray[i] = false;
    NextActiveArray[i] = false;
  }
  for(long i = 0; i < batch_size; i++) {
    ShortestPathLen[(IdxType)batch_size * (IdxType)vecQueries[i] + (IdxType)i] = 0;
  }
  for(size_t i = 0; i < batch_size; i++) {
    CurrActiveArray[(IdxType)vecQueries[i] * (IdxType)batch_size + (IdxType)i] = true;
  }

  vertexSubset Frontier(n, frontier);

  // for profiling
  long iteration = 0;
  size_t totalActivated = 0;
  size_t totalNoOverlap = 0;
  vector<pair<size_t, double>> affinity_tracking;
  vector<pair<size_t, double>> affinity_tracking_one;
  vector<pair<size_t, double>> affinity_tracking_only;
  size_t* overlaps = pbbs::new_array<size_t>(batch_size);

  for (int i = 0; i < batch_size; i++) {
    overlaps[i] = 0;
  }

  vector<double> overlap_scores;
  size_t accumulated_overlap = 0;
  size_t accumulated_overlap_one = 0;
  size_t accumulated_overlap_only= 0;
  size_t peak_activation = 0;
  int peak_iter = 0;

  // vector<long> frontier_iterations;
  // vector<long> overlapped_iterations;
  // vector<long> accumulated_overlapped_iterations;
  // vector<long> total_activated_iterations;

  while(!Frontier.isEmpty()){
    iteration++;
    totalActivated += Frontier.size();
    // cout << Frontier.size() << endl;
    // profiling
    if (should_profile) {
      
    }

    // mode: no_dense, remove_duplicates (for batch size > 1)
    vertexSubset output = edgeMap(G, Frontier, DJ_F(ShortestPathLen, CurrActiveArray, NextActiveArray, batch_size), -1, no_dense|remove_duplicates);

    Frontier.del();
    Frontier = output;

    Frontier.toDense();
    bool* new_d = Frontier.d;
    Frontier.d = nullptr;
    vertexSubset Frontier_new(n, new_d);
    Frontier.del();
    Frontier = Frontier_new;

    std::swap(CurrActiveArray, NextActiveArray);
    parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
      NextActiveArray[i] = false;
    }
  }
  cout << "Total iterations: " << iteration << endl;
  // cout << endl;
  // profiling
  if (should_profile) {
    
  }

#ifdef OUTPUT 
  for (int i = 0; i < batch_size; i++) {
    long start = vecQueries[i];
    char outFileName[300];
    sprintf(outFileName, "SSSP_base_output_src%ld.%ld.%ld.out", start, edge_count, batch_size);
    FILE *fp;
    fp = fopen(outFileName, "w");
    for (long j = 0; j < n; j++)
      fprintf(fp, "%ld %d\n", j, ShortestPathLen[j * batch_size + i]);
    fclose(fp);
  }
#endif

  Frontier.del();
  pbbs::delete_array(ShortestPathLen, totalNumVertices);
  pbbs::delete_array(CurrActiveArray, totalNumVertices);
  pbbs::delete_array(NextActiveArray, totalNumVertices);
  pbbs::delete_array(overlaps, batch_size);
  return make_pair(totalActivated, totalNoOverlap);
}

template <class vertex>
pair<size_t, size_t> Compute_Separate(graph<vertex>& G, std::vector<long> vecQueries, commandLine P, bool should_profile) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();

  cout << "initializing\n";
  intE** vals = pbbs::new_array<intE*>(batch_size);
  bool** frontiers = pbbs::new_array<bool*>(n);
  for (int i = 0; i < batch_size; i++) {
    vals[i] = pbbs::new_array<intE>(n);
    frontiers[i] = pbbs::new_array<bool>(n);
    parallel_for(IdxType j = 0; j < n; j++) {
      vals[i][j] = (intE)MAXPATH;
      frontiers[i][j] = false;
    }
    vals[i][(IdxType)vecQueries[i]] = 0;
    frontiers[i][(IdxType)vecQueries[i]] = true;
  }

  // for profiling
  long iteration = 0;
  size_t totalActivated = 0;
  size_t totalNoOverlap = 0;

  bool isConverged = true;
  vertexSubset** vecFs = new vertexSubset*[batch_size];
  // vector<vertexSubset*> vecFs(batch_size);
  for (int i = 0; i < batch_size; i++) {
    // vertexSubset Frontier_new(n, frontiers[i]);

    vecFs[i] = new vertexSubset(n, frontiers[i]);
    // vecFs.push_back(&Frontier_new);
    isConverged = isConverged && vecFs[i]->isEmpty();
  }
  cout << "finished initializing\n";

  while (!isConverged) {
    iteration++;
    // cout << "iteration: " << iteration << endl;
    parallel_for(int j = 0; j < batch_size; j++) {
      vertexSubset output = edgeMap(G, *vecFs[j], DJ_SINGLE_F(vals[j]), -1, no_dense|remove_duplicates);
      vecFs[j]->del();
      *vecFs[j] = output;
    }
    bool isEmpty = true;
    for (int i = 0; i < batch_size; i++) {
      isEmpty = isEmpty && vecFs[i]->isEmpty();
    }
    isConverged = isEmpty;
  }
  cout << "Total iterations: " << iteration << endl;

#ifdef OUTPUT 
  for (int i = 0; i < batch_size; i++) {
    long start = vecQueries[i];
    char outFileName[300];
    sprintf(outFileName, "SSSP_separate_output_src%ld.%ld.%ld.out", start, edge_count, batch_size);
    FILE *fp;
    fp = fopen(outFileName, "w");
    for (long j = 0; j < n; j++)
      fprintf(fp, "%ld %d\n", j, vals[i][j]);
    fclose(fp);
  }
#endif

  // Frontier.del();
  for (int i = 0; i < batch_size; i++) {
    vecFs[i]->del();
    pbbs::delete_array(vals[i], n);
  }

  return make_pair(totalActivated, totalNoOverlap);
}

template <class vertex>
pair<size_t, size_t> Compute_Base_Dynamic(graph<vertex>& G, std::vector<long> vecQueries, queue<long>& queryQueue, commandLine P, bool should_profile) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size;
  intE* ShortestPathLen = pbbs::new_array<intE>(totalNumVertices);
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
    ShortestPathLen[i] = (intE)MAXPATH;
    CurrActiveArray[i] = false;
    NextActiveArray[i] = false;
  }
  for(long i = 0; i < batch_size; i++) {
    ShortestPathLen[(IdxType)batch_size * (IdxType)vecQueries[i] + (IdxType)i] = 0;
  }
  for(size_t i = 0; i < batch_size; i++) {
    CurrActiveArray[(IdxType)vecQueries[i] * (IdxType)batch_size + (IdxType)i] = true;
  }

  vertexSubset Frontier(n, frontier);

  
  long iteration = 0;
  size_t totalActivated = 0;
  size_t totalNoOverlap = 0;

  // for async batching
  int slots_parameter = P.getOptionIntValue("-slots", 1);
  int empty_slots = 0;
  vector<long> next_batch;
  vector<long> lastQueries;
  bool _flag_check = true;
  bool* finish_status = pbbs::new_array<bool>(batch_size);
  for(long i = 0; i < batch_size; i++) {
    finish_status[i] = true;
  }
  
  vector<long> queuedQueries;
  long queued_index = 0;
  while (!queryQueue.empty()) {
    queuedQueries.push_back(queryQueue.front());
    queryQueue.pop();
  }

  while(!Frontier.isEmpty() || queued_index != queuedQueries.size()){
    iteration++;
    // cout << "iteration: " << iteration << endl;
    totalActivated += Frontier.size();

    // mode: no_dense, remove_duplicates (for batch size > 1)
    vertexSubset output = edgeMap(G, Frontier, DJ_F(ShortestPathLen, CurrActiveArray, NextActiveArray, batch_size), -1, no_dense|remove_duplicates);

    Frontier.del();
    Frontier = output;

    std::swap(CurrActiveArray, NextActiveArray);
    parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
      NextActiveArray[i] = false;
    }

    bool* is_batch_finished = pbbs::new_array<bool>(batch_size);
    for(size_t i = 0; i < batch_size; i++) {
      is_batch_finished[i] = false;
    }
    empty_slots = 0;
    if (queued_index != queuedQueries.size()) { // there are more queries to process
      vector<long> empty_slots_vec;
      for (int i = 0; i < batch_size; i++) {
        bool is_finished = true;
        long active_cnt = 0;
        parallel_for(size_t index = 0; index < n; index++) {
          if (CurrActiveArray[index * batch_size + i]) {
            // is_finished = false;
            pbbs::fetch_and_add(&active_cnt, 1);
          }
        }
        if (active_cnt == 0) {
          empty_slots++;
          empty_slots_vec.push_back(i);
          if (finish_status[i]) {
            finish_status[i] = false;
          }
        }
      }
      if (empty_slots == batch_size) {
        cout << "batch is empty, iteration " << iteration << endl;
      }

      if (empty_slots >= slots_parameter || (queuedQueries.size()-queued_index <= empty_slots)) {
        for (int i = 0; i < empty_slots_vec.size(); i++) {
          if (queued_index < queuedQueries.size()) {
            long nextQ = queuedQueries[queued_index];
            queued_index++;
            int index_in_batch = empty_slots_vec[i];
            #ifdef OUTPUT 
              long start = vecQueries[index_in_batch];
              char outFileName[300];
              sprintf(outFileName, "SSSP_base_async_output_src%ld.%ld.%ld.out", start, edge_count, batch_size);
              FILE *fp;
              fp = fopen(outFileName, "w");
              for (long j = 0; j < n; j++)
                fprintf(fp, "%ld %d\n", j, ShortestPathLen[j * batch_size + index_in_batch]);
              fclose(fp);
            #endif
            vecQueries[index_in_batch] = nextQ;
            finish_status[index_in_batch] = true;

            CurrActiveArray[nextQ * batch_size + index_in_batch] = true;
            parallel_for(size_t index = 0; index < n; index++) {
              ShortestPathLen[index * batch_size + index_in_batch] = (intE)MAXPATH;
            }
            ShortestPathLen[nextQ * batch_size + index_in_batch] = 0;

            Frontier.toDense();
            bool* new_d = pbbs::new_array<bool>(n);
            parallel_for(size_t ii = 0; ii < n; ii++) {
              new_d[ii] = Frontier.d[ii];
            }
            new_d[nextQ] = true;
            vertexSubset Frontier_new(n, new_d);
            Frontier.del();
            Frontier = Frontier_new;

            // cout << "iteration: " << iteration << ": inserting " << nextQ << " into the batch\n";
            // cout << "current empty slots: " << empty_slots << endl;
          }
        }
      }
    }
    // if (queued_index >= queuedQueries.size() && lastQueries.size() == 0) {
    //   for (int i = 0; i < batch_size; i++) {

    //   }
    // }
  }

  // profiling
  if (should_profile) {
    
  }

#ifdef OUTPUT 
  for (int i = 0; i < batch_size; i++) {
    long start = vecQueries[i];
    char outFileName[300];
    sprintf(outFileName, "SSSP_base_async_output_src%ld.%ld.%ld.out", start, edge_count, batch_size);
    FILE *fp;
    fp = fopen(outFileName, "w");
    for (long j = 0; j < n; j++)
      fprintf(fp, "%ld %d\n", j, ShortestPathLen[j * batch_size + i]);
    fclose(fp);
  }
#endif

  Frontier.del();
  pbbs::delete_array(ShortestPathLen, totalNumVertices);
  pbbs::delete_array(CurrActiveArray, totalNumVertices);
  pbbs::delete_array(NextActiveArray, totalNumVertices);
  return make_pair(totalActivated, totalNoOverlap);
}

template <class vertex>
pair<size_t, size_t> Compute_Base_Skipping(graph<vertex>& G, std::vector<long> vecQueries, commandLine P, int skipIter, bool should_profile) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size;
  intE* ShortestPathLen = pbbs::new_array<intE>(totalNumVertices);
  bool* frontier = pbbs::new_array<bool>(n);
  parallel_for(size_t i = 0; i < n; i++) {
    frontier[i] = false;
  }
  for(long i = 0; i < batch_size; i++) {
    frontier[vecQueries[i]] = true;
  }
  parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
    ShortestPathLen[i] = (intE)MAXPATH;
  }
  for(long i = 0; i < batch_size; i++) {
    ShortestPathLen[(IdxType)batch_size * (IdxType)vecQueries[i] + (IdxType)i] = 0;
  }

  vertexSubset Frontier(n, frontier);

  // for profiling
  long iteration = 0;
  size_t totalActivated = 0;
  size_t totalNoOverlap = 0;
  size_t total_edges = 0;
  while(!Frontier.isEmpty()){
    iteration++;
    totalActivated += Frontier.size();
    cout << "iteration: " << Frontier.size() << ":: " << totalActivated << endl;
    // profile edge activations.
    Frontier.toDense();
    for (IdxType i = 0; i < n; i++) {
      if (Frontier.d[i]) {
        long temp_degree =  G.V[i].getOutDegree();
        total_edges += temp_degree;
      }
    }

    // mode: no_dense, remove_duplicates (for batch size > 1)
    if (iteration > skipIter) {
      vertexSubset output = edgeMap(G, Frontier, DJ_SKIP_F(ShortestPathLen, batch_size), -1, no_dense|remove_duplicates);
      Frontier.del();
      Frontier = output;

      Frontier.toDense();
      bool* new_d = Frontier.d;
      Frontier.d = nullptr;
      vertexSubset Frontier_new(n, new_d);
      Frontier.del();
      Frontier = Frontier_new;
    }
  }
  cout << "Total iterations: " << iteration << endl;


#ifdef OUTPUT 
  for (int i = 0; i < batch_size; i++) {
    long start = vecQueries[i];
    char outFileName[300];
    sprintf(outFileName, "SSSP_base_skipping_output_src%ld.%ld.%ld.out", start, edge_count, batch_size);
    FILE *fp;
    fp = fopen(outFileName, "w");
    for (long j = 0; j < n; j++)
      fprintf(fp, "%ld %d\n", j, ShortestPathLen[j * batch_size + i]);
    fclose(fp);
  }
#endif

  Frontier.del();
  pbbs::delete_array(ShortestPathLen, totalNumVertices);
  return make_pair(total_edges, totalNoOverlap);
}

template <class vertex>
pair<size_t, size_t> Compute_Delay(graph<vertex>& G, std::vector<long> vecQueries, commandLine P, std::vector<int> defer_vec, bool should_profile) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size;
  intE* ShortestPathLen = pbbs::new_array<intE>(totalNumVertices);
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
  // for(long i = 0; i < batch_size; i++) {
  //   frontier[vecQueries[i]] = true;
  // }
  parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
    ShortestPathLen[i] = (intE)MAXPATH;
    CurrActiveArray[i] = false;
    NextActiveArray[i] = false;
  }
  for(long i = 0; i < batch_size; i++) {
    ShortestPathLen[(IdxType)batch_size * (IdxType)vecQueries[i] + (IdxType)i] = 0;
  }
  for(size_t i = 0; i < batch_size; i++) {
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
    vertexSubset output = edgeMap(G, Frontier, DJ_F(ShortestPathLen, CurrActiveArray, NextActiveArray, batch_size), -1, no_dense|remove_duplicates);

    Frontier.del();
    Frontier = output;

    // Checking for delayed queries.
    // TODO: possibly optimization.
    // timer t_checking;
    // t_checking.start();
    // Frontier.toDense();
    // bool* new_d = pbbs::new_array<bool>(n);
    // parallel_for(size_t i = 0; i < n; i++) {
    //   new_d[i] = Frontier.d[i];
    // }

    
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
    // t_checking.stop();
    // t_checking.reportTotal("checking delaying");

    std::swap(CurrActiveArray, NextActiveArray);
    parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
      NextActiveArray[i] = false;
    }
  }

#ifdef OUTPUT 
  for (int i = 0; i < batch_size; i++) {
    long start = vecQueries[i];
    char outFileName[300];
    sprintf(outFileName, "SSSP_delay_output_src%ld.%ld.%ld.out", start, edge_count, batch_size);
    FILE *fp;
    fp = fopen(outFileName, "w");
    for (long j = 0; j < n; j++)
      fprintf(fp, "%ld %d\n", j, ShortestPathLen[j * batch_size + i]);
    fclose(fp);
  }
#endif

  Frontier.del();
  pbbs::delete_array(CurrActiveArray, totalNumVertices);
  pbbs::delete_array(NextActiveArray, totalNumVertices);
  pbbs::delete_array(ShortestPathLen, totalNumVertices);
  return make_pair(totalActivated, 0);
}

template <class vertex>
pair<size_t, size_t> Compute_Delay_Skipping(graph<vertex>& G, std::vector<long> vecQueries, commandLine P, std::vector<int> defer_vec, bool should_profile) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size;
  intE* ShortestPathLen = pbbs::new_array<intE>(totalNumVertices);
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
  // for(long i = 0; i < batch_size; i++) {
  //   frontier[vecQueries[i]] = true;
  // }
  parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
    ShortestPathLen[i] = (intE)MAXPATH;
  }
  for(long i = 0; i < batch_size; i++) {
    ShortestPathLen[(IdxType)batch_size * (IdxType)vecQueries[i] + (IdxType)i] = 0;
  }

  vertexSubset Frontier(n, frontier);

  // for profiling
  long iteration = 0;
  size_t totalActivated = 0;
  size_t total_edges = 0;
  while(!Frontier.isEmpty()){
    iteration++;
    totalActivated += Frontier.size();
    cout << "iteration: " << Frontier.size() << ":: " << totalActivated << endl;
    // profile edge activations.
    Frontier.toDense();
    for (IdxType i = 0; i < n; i++) {
      if (Frontier.d[i]) {
        long temp_degree =  G.V[i].getOutDegree();
        total_edges += temp_degree;
      }
    }

    // mode: no_dense, remove_duplicates (for batch size > 1)
    vertexSubset output = edgeMap(G, Frontier, DJ_SKIP_F(ShortestPathLen, batch_size), -1, no_dense|remove_duplicates);
    Frontier.del();
    Frontier = output;
    
    Frontier.toDense();
    bool* new_d = Frontier.d;
    Frontier.d = nullptr;
    for(long i = 0; i < batch_size; i++) {
      if (defer_vec[i] == iteration) {
        new_d[vecQueries[i]] = true;
        // NextActiveArray[(IdxType)vecQueries[i] * (IdxType)batch_size + (IdxType)i] = true;
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
    sprintf(outFileName, "SSSP_delay_output_src%ld.%ld.%ld.out", start, edge_count, batch_size);
    FILE *fp;
    fp = fopen(outFileName, "w");
    for (long j = 0; j < n; j++)
      fprintf(fp, "%ld %d\n", j, ShortestPathLen[j * batch_size + i]);
    fclose(fp);
  }
#endif

  Frontier.del();
  pbbs::delete_array(ShortestPathLen, totalNumVertices);
  return make_pair(total_edges, 0);
}