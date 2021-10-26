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

struct DJ_Chunk_F {
  intE* ShortestPathLen;
  long BatchSize;

  DJ_Chunk_F(intE* _ShortestPathLen, long _BatchSize) : 
    ShortestPathLen(_ShortestPathLen), BatchSize(_BatchSize) {}
  
  inline bool update (uintE s, uintE d, intE edgeLen) { //Update
    bool ret = false;
    IdxType s_begin = s * BatchSize;
    IdxType d_begin = d * BatchSize;

    for (long j = 0; j < BatchSize; j++) {
      intE newValue = ShortestPathLen[s_begin + j] + edgeLen;
      if (ShortestPathLen[d_begin + j] > newValue) {
        // Visited[d_begin + j] = true;
        ShortestPathLen[d_begin + j] = newValue;
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
        intE newValue = ShortestPathLen[s_begin + j] + edgeLen;
        if (writeMin(&ShortestPathLen[d_begin + j], newValue)) {
          ret = true;
        }
    }
    return ret;
  }
  //cond function checks if vertex has been visited yet
  inline bool cond (uintE d) { return cond_true(d); } 
};

struct DJ_SKIP_F {
  intE* ShortestPathLen;
  bool* CurrActiveArray;
  bool* NextActiveArray;
  long BatchSize;

  DJ_SKIP_F(intE* _ShortestPathLen, bool* _CurrActiveArray, bool* _NextActiveArray, long _BatchSize) : 
    ShortestPathLen(_ShortestPathLen), CurrActiveArray(_CurrActiveArray), NextActiveArray(_NextActiveArray), BatchSize(_BatchSize) {}
  
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

  while(!Frontier.isEmpty()){
    iteration++;
    totalActivated += Frontier.size();

    // mode: no_dense, remove_duplicates (for batch size > 1)
    vertexSubset output = edgeMap(G, Frontier, DJ_F(ShortestPathLen, CurrActiveArray, NextActiveArray, batch_size), -1, no_dense|remove_duplicates);

    Frontier.del();
    Frontier = output;

    std::swap(CurrActiveArray, NextActiveArray);
    parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
      NextActiveArray[i] = false;
    }
  }
  // cout << endl;

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
pair<size_t, size_t> Compute_Chunk(graph<vertex>& G, std::vector<long> vecQueries, const vector<set<long>>& Tables, const long* vtx2chunk, commandLine P, bool should_profile) {
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

  // vector<set<long>> Tables;
  // Tables.insert(Tables.end(), C_Set.begin(), C_Set.end());
  long chunk_size = Tables.size();
  // long* vtx2chunk = pbbs::new_array<long>(n);
  // // long* vtx2chunk = newA(long,n);
  // parallel_for(long i = 0; i < n; i++) {
  //   vtx2chunk[i] = -1;
  // }
  // parallel_for(long i = 0; i < chunk_size; i++) {
  //   set<long> tmp_chunk = Tables[i];
  //   // cout << "Chunk " << i << endl;
  //   for (auto e : tmp_chunk) {
  //     // cout << "\t " << e << endl;
  //     vtx2chunk[e] = i;
  //   }
  // }
  cout << "chunk_size: " << chunk_size << endl;
  long* chunkCnt = pbbs::new_array<long>(chunk_size);
  parallel_for(long i = 0; i < chunk_size; i++) {
    chunkCnt[i] = 0;
  }
  cout << chunkCnt << endl;
  // vector<long> chunkCnt(chunk_size);
  cout << "Query: " << vecQueries[0] << endl;
  while(!Frontier.isEmpty()){
    iteration++;
    cout << "iteration: " << iteration << endl;
    totalActivated += Frontier.size();
    // cout << "Frontier size: " << Frontier.size() << endl;
    // should find the chunk to be processed.
    // Such a chunk is needed by most number of jobs. 
    Frontier.toDense();
    bool* f_old = Frontier.d;
    Frontier.d = nullptr;
    Frontier.del();

    // bool* f_new = pbbs::new_array<bool>(n);
    // bool* f_new = newA(bool,n);
    // bool* f_next = pbbs::new_array<bool>(n);
    bool* f_next = newA(bool,n);
    // bool* f_remains = pbbs::new_array<bool>(n);
    bool* f_remains = newA(bool,n);
    // long* chunkCnt = pbbs::new_array<long>(chunk_size);
    // long* chunkCnt = newA(long,chunk_size);
    // vector<long> chunkCnt(chunk_size);
    // parallel_for(long i = 0; i < chunk_size; i++) {
    //   chunkCnt[i] = 0;
    // }
    {parallel_for(long i = 0; i < n; i++) {
      // f_new[i] = f_old[i];
      f_next[i] = false;
      f_remains[i] = false;
      if (f_old[i]) {
        // if (vtx2chunk[i] == -1) cout << "-1, the error!!\n";
        writeAdd(&chunkCnt[vtx2chunk[i]],1l);
      }
    }}
    // cout << "Chunk " << distance(chunkCnt, max_element(chunkCnt, chunkCnt + chunk_size))  << " is required most." << endl;
    // long load_chunk_id = distance(chunkCnt, max_element(chunkCnt, chunkCnt + chunk_size));
    // cout << "index = " << load_chunk_id << endl;
    long* max_chunk_needed = max_element(chunkCnt, chunkCnt + chunk_size);
    long load_chunk_id = max_chunk_needed-chunkCnt;
    // cout << "max count needed " << *max_chunk_needed << ", index = " << max_chunk_needed-chunkCnt << endl;
    // long load_chunk_id = std::distance(chunkCnt.begin(),std::max_element(chunkCnt.begin(), chunkCnt.end()));

    // cout << "\tDEBUG pos a-1: " << endl;
    {parallel_for(long i = 0; i < n; i++) {
      if (f_old[i]) {
        // if (vtx2chunk[i] == -1) cout << "-1, the error!!\n";
        if (vtx2chunk[i] != load_chunk_id) {
          f_old[i] = false;
          f_remains[i] = true;
        }
      }
    }}
    // cout << "\tDEBUG pos a-2: " << endl;
    // Frontier.d = nullptr;
    vertexSubset Frontier_new(n, f_old);
    // cout << "\tDEBUG pos a-3: " << endl;
    // 
    // cout << "start processing chunks... iteration " << iteration << "\n";
    while (!Frontier_new.isEmpty()) {
      // Frontier_new.toDense();
      // for (long i = 0; i < n;i++) {
      //   if (Frontier_new.d[i]) {
      //     cout << i << endl;
      //   }
      // }
      // cout << "\tFrontier_new size: " << Frontier_new.size() << endl;
      vertexSubset output = edgeMap(G, Frontier_new, DJ_Chunk_F(ShortestPathLen, batch_size), -1, no_dense|remove_duplicates);
      // cout << "\t\t DEBUG pos 1: " << endl;
      // bool* f_next_chunk = newA(bool,n);
      // cout << "\t\t output size: " << output.size() << endl;
      output.toDense();
      // cout << "\t\t DEBUG pos 1-2: " << endl;
      bool* new_output = output.d;
      {parallel_for(long i = 0; i < n; i++) {
        // f_next_chunk[i] = output.d[i];
        if (new_output[i]) {
          // cout << "new_output[" << i << "]\n";
          // if (vtx2chunk[i] == -1) cout << "-1, the error!!\n";
          if (vtx2chunk[i] != load_chunk_id) {
            // cout << "load_chunk_id " << load_chunk_id << "\n";
            f_next[i] = true;
            new_output[i] = false;
          }
        }
      }}
      // cout << "\t\t DEBUG pos 2: " << endl;
      
      // cout << "\t\t DEBUG pos 3: " << endl;
      // vertexSubset Next_Chunk(n, f_next_chunk);
      // cout << "\t\t DEBUG pos 4: " << endl;
      vertexSubset Next_Chunk(n, new_output);
      output.d = nullptr;
      output.del();
      // cout << "\tNext_Chunk size: " << Next_Chunk.size() << endl;
      if (Next_Chunk.isEmpty())
      {
        // cout << "\t\t\tNext_Chunk is empty!\n";
        // current chunk is finished, should proceed next (remaining) chunk
        // bool* f_new_chunk = pbbs::new_array<bool>(n);
        // long* _chunkCnt = pbbs::new_array<long>(chunk_size);
        // bool* f_new_chunk = newA(bool,n);
        Frontier_new.toDense();
        bool* f_new_chunk = Frontier_new.d;
        Frontier_new.d = nullptr;
        // long* _chunkCnt = newA(long,chunk_size);
        // parallel_for(long i = 0; i < chunk_size; i++) {
        //   _chunkCnt[i] = 0;
        // }
        parallel_for(long i = 0; i < chunk_size; i++) {
          chunkCnt[i] = 0;
        }
        // cout << "\t\t DEBUG pos 5: " << endl;
        {parallel_for(long i = 0; i < n; i++) {
          // f_new_chunk[i] = f_remains[i];
          if (f_remains[i]) {
            writeAdd(&chunkCnt[vtx2chunk[i]],1l);
          }
        }}
        // cout << "\t\t DEBUG pos 6: " << endl;
        // long _load_chunk_id = distance(chunkCnt, max_element(chunkCnt, chunkCnt + chunk_size));
        // cout << "index = " << _load_chunk_id << endl;
        long* _max_chunk_needed = max_element(chunkCnt, chunkCnt + chunk_size);
        long _load_chunk_id = _max_chunk_needed-chunkCnt;
        // cout << "max count needed " << *_max_chunk_needed << ", index = " << _max_chunk_needed-chunkCnt << endl;
        // long _load_chunk_id = std::distance(chunkCnt.begin(),std::max_element(chunkCnt.begin(), chunkCnt.end()));
        {parallel_for(long i = 0; i < n; i++) {
          f_new_chunk[i] = false;
          if (f_remains[i] && (vtx2chunk[i] == _load_chunk_id)) {
            f_new_chunk[i] = true;
            f_remains[i] = false;
          }
        }}
        // cout << "\t\t DEBUG pos 6-2: " << endl;
        // pbbs::delete_array(_chunkCnt, chunk_size);
        // free(_chunkCnt);
        // cout << "\t\t DEBUG pos 7: " << endl;
        vertexSubset _Next_Chunk(n, f_new_chunk);
        // cout << "\t\t DEBUG pos 8: " << endl;
        
        Frontier_new.del();
        // cout << "\t\t DEBUG pos 9: " << endl;
        Frontier_new = _Next_Chunk;
        // cout << "\t\t DEBUG pos 10: " << endl;
        Next_Chunk.del();
        // cout << "\t\tFinish Next_Chunk is empty!\n";
      } else {
        // cout << "\t\tNext_Chunk is NOT empty!\n";
        Frontier_new.del();
        Frontier_new = Next_Chunk;
      }
      // cout << "\t\t DEBUG pos 11: " << endl;
    }
    // cout << "\t\t DEBUG pos 12: " << endl;
    // bool* f_new_d = pbbs::new_array<bool>(n);
    vertexSubset Next_Frontier(n, f_next);
    // cout << "\t\t DEBUG pos 13: " << endl;
    
    // cout << "\t\t DEBUG pos 14: " << endl;

    Frontier = Next_Frontier;
    // cout << "\t\t DEBUG pos 15: " << endl;
    // pbbs::delete_array(chunkCnt, chunk_size);
    // free(chunkCnt);
    // cout << "== finish one iteration\n";
    parallel_for(long i = 0; i < chunk_size; i++) {
      chunkCnt[i] = 0;
    }
  }
  cout << "========finished query(ies)========" << endl;

#ifdef OUTPUT 
  for (int i = 0; i < batch_size; i++) {
    long start = vecQueries[i];
    char outFileName[300];
    sprintf(outFileName, "SSSP_chunk_output_src%ld.%ld.%ld.out", start, edge_count, batch_size);
    FILE *fp;
    fp = fopen(outFileName, "w");
    for (long j = 0; j < n; j++)
      fprintf(fp, "%ld %d\n", j, ShortestPathLen[j * batch_size + i]);
    fclose(fp);
  }
#endif
cout << "========cleaning========" << endl;
  Frontier.del();
  // free(vtx2chunk);
  // pbbs::delete_array(vtx2chunk, n);
  cout << chunkCnt << endl;
  pbbs::delete_array(chunkCnt, chunk_size);
  pbbs::delete_array(ShortestPathLen, totalNumVertices);
cout << "========finished cleaning========" << endl;
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

  while(!Frontier.isEmpty()){
    iteration++;
    totalActivated += Frontier.size();

    // mode: no_dense, remove_duplicates (for batch size > 1)
    if (iteration > skipIter) {
      vertexSubset output = edgeMap(G, Frontier, DJ_SKIP_F(ShortestPathLen, CurrActiveArray, NextActiveArray, batch_size), -1, no_dense|remove_duplicates);
      Frontier.del();
      Frontier = output;

      Frontier.toDense();
      bool* new_d = Frontier.d;
      Frontier.d = nullptr;
      vertexSubset Frontier_new(n, new_d);
      Frontier.del();
      Frontier = Frontier_new;
    } else {
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
  }


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
  pbbs::delete_array(CurrActiveArray, totalNumVertices);
  pbbs::delete_array(NextActiveArray, totalNumVertices);
  return make_pair(totalActivated, totalNoOverlap);
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
    vertexSubset output = edgeMap(G, Frontier, DJ_SKIP_F(ShortestPathLen, CurrActiveArray, NextActiveArray, batch_size), -1, no_dense|remove_duplicates);
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
  pbbs::delete_array(CurrActiveArray, totalNumVertices);
  pbbs::delete_array(NextActiveArray, totalNumVertices);
  pbbs::delete_array(ShortestPathLen, totalNumVertices);
  return make_pair(totalActivated, 0);
}