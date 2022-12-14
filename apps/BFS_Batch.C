// This code is part of the project "Glign: Taming Misaligned Graph Traversals in Concurrent Graph Processing"
// Copyright (c) 2022 Xizhe Yin
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

struct BFSLV_SKIP_F {
  uintE* Levels;
  long BatchSize;

  BFSLV_SKIP_F(uintE* _Levels, long _BatchSize) : 
    Levels(_Levels), BatchSize(_BatchSize) {}
  
  inline bool update (uintE s, uintE d, intE edgeLen=1) { //Update
    bool ret = false;
    IdxType s_begin = s * BatchSize;
    IdxType d_begin = d * BatchSize;

    for (long j = 0; j < BatchSize; j++) {
      uintE newValue = Levels[s_begin + j] + 1;
      if (Levels[d_begin + j] > newValue) {
        Levels[d_begin + j] = newValue;
        ret = true;
      }
    }
    return ret;
  }
  
  inline bool updateAtomic (uintE s, uintE d, intE edgeLen=1){ //atomic version of Update
    bool ret = false;
    IdxType s_begin = s * BatchSize;
    IdxType d_begin = d * BatchSize;
    for (long j = 0; j < BatchSize; j++) {
      uintE newValue = Levels[s_begin + j] + 1;
      if (writeMin(&Levels[d_begin + j], newValue)) {
        ret = true;
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

  // for profiling
  long iteration = 0;
  size_t totalActivated = 0;

  while(!Frontier.isEmpty()){
    iteration++;
    totalActivated += Frontier.size();
    // cout << "iteration: " << Frontier.size() << ": " << totalActivated << endl;
    // mode: no_dense, remove_duplicates (for batch size > 1)
    vertexSubset output = edgeMap(G, Frontier, BFSLV_F(Levels, CurrActiveArray, NextActiveArray, batch_size), -1, no_dense|remove_duplicates);

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
    sprintf(outFileName, "BFSLV_base_output_src%ld.%ld.%ld.out", start, edge_count, batch_size);
    FILE *fp;
    fp = fopen(outFileName, "w");
    for (long j = 0; j < n; j++)
      fprintf(fp, "%ld %d\n", j, Levels[j * batch_size + i]);
    fclose(fp);
  }
#endif

  Frontier.del();
  pbbs::delete_array(CurrActiveArray, totalNumVertices);
  pbbs::delete_array(NextActiveArray, totalNumVertices);
  pbbs::delete_array(Levels, totalNumVertices);
  return make_pair(totalActivated, 0);
}

template <class vertex>
pair<size_t, size_t> Compute_Base_Dynamic(graph<vertex>& G, std::vector<long> vecQueries, queue<long>& queryQueue, commandLine P, bool should_profile) {
  return make_pair(0,0);
}

template <class vertex>
pair<size_t, size_t> Compute_Separate(graph<vertex>& G, std::vector<long> vecQueries, commandLine P, bool should_profile) {
  return make_pair(0,0);
}

template <class vertex>
pair<size_t, size_t> Compute_Base_Skipping(graph<vertex>& G, std::vector<long> vecQueries, commandLine P, bool should_profile) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size;
  uintE* Levels = pbbs::new_array<uintE>(totalNumVertices);
  bool* frontier = pbbs::new_array<bool>(n);
  parallel_for(size_t i = 0; i < n; i++) {
    frontier[i] = false;
  }
  for(long i = 0; i < batch_size; i++) {
    frontier[vecQueries[i]] = true;
  }
  parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
    Levels[i] = (uintE)MAXLEVEL;
  }
  for(long i = 0; i < batch_size; i++) {
    Levels[(IdxType)batch_size * (IdxType)vecQueries[i] + (IdxType)i] = 0;
  }

  vertexSubset Frontier(n, frontier);

  // for profiling
  long iteration = 0;
  size_t totalActivated = 0;

  while(!Frontier.isEmpty()){
    iteration++;
    totalActivated += Frontier.size();
    // cout << "iteration: " << iteration << ": " << Frontier.size() << endl;
    // mode: no_dense, remove_duplicates (for batch size > 1)
    vertexSubset output = edgeMap(G, Frontier, BFSLV_SKIP_F(Levels, batch_size), -1, no_dense|remove_duplicates);
    
    Frontier.del();
    Frontier = output;

    Frontier.toDense();
    bool* new_d = Frontier.d;
    Frontier.d = nullptr;
    vertexSubset Frontier_new(n, new_d);
    Frontier.del();
    Frontier = Frontier_new;
  }

#ifdef OUTPUT 
  for (int i = 0; i < batch_size; i++) {
    long start = vecQueries[i];
    char outFileName[300];
    sprintf(outFileName, "BFSLV_base_output_src%ld.%ld.%ld.out", start, edge_count, batch_size);
    FILE *fp;
    fp = fopen(outFileName, "w");
    for (long j = 0; j < n; j++)
      fprintf(fp, "%ld %d\n", j, Levels[j * batch_size + i]);
    fclose(fp);
  }
#endif

  Frontier.del();
  pbbs::delete_array(Levels, totalNumVertices);
  return make_pair(totalActivated, 0);
}

template <class vertex>
pair<size_t, size_t> Compute_Delay(graph<vertex>& G, std::vector<long> vecQueries, commandLine P, std::vector<int> defer_vec, bool should_profile) {
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
  // for delaying initialization
  for(long i = 0; i < batch_size; i++) {
    if (defer_vec[i] == 0 || defer_vec[i] > 1000) {
      frontier[vecQueries[i]] = true;
    }
  }
  // for(long i = 0; i < batch_size; i++) {
  //   frontier[vecQueries[i]] = true;
  // }
  parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
    Levels[i] = (uintE)MAXLEVEL;
    CurrActiveArray[i] = false;
    NextActiveArray[i] = false;
  }
  for(long i = 0; i < batch_size; i++) {
    Levels[(IdxType)batch_size * (IdxType)vecQueries[i] + (IdxType)i] = 0;
  }
  parallel_for(size_t i = 0; i < batch_size; i++) {
    if (defer_vec[i] == 0 || defer_vec[i] > 1000) {
      CurrActiveArray[(IdxType)vecQueries[i] * (IdxType)batch_size + (IdxType)i] = true;
    }
  }

  vertexSubset Frontier(n, frontier);

  // for profiling
  long iteration = 0;
  size_t totalActivated = 0;
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

    // profiling
    if (should_profile) {
      // if (Frontier.size() > peak_activation) {
      //   peak_activation = Frontier.size();
      //   peak_iter = iteration;
      // }
      // bool* overlap_set = pbbs::new_array<bool>(n); // activated for all queries
      // bool* overlap_set_one = pbbs::new_array<bool>(n); // only activated for at least half of queries
      // bool* overlap_set_only = pbbs::new_array<bool>(n);
      // parallel_for(size_t i = 0; i < n; i++) {
      //   overlap_set[i] = true;
      //   overlap_set_one[i] = false;
      //   overlap_set_only[i] = false;
      // }
      // parallel_for(size_t index = 0; index < n; index++) {
      //   int tmp_flag = 0;
      //   for (int i = 0; i < batch_size; i++) {
      //     // if the vertex is activated for all queries...
      //     overlap_set[index] = overlap_set[index] && CurrActiveArray[index * batch_size + i];
      //     if (CurrActiveArray[index * batch_size + i]) {
      //       tmp_flag++;
      //     }
      //   }
      //   // activated for at least half of the queries.
      //   if (tmp_flag >= batch_size/2) {
      //     overlap_set_one[index] = true;
      //   }
      //   // only one is activated.
      //   if (tmp_flag == 1) {
      //     overlap_set_only[index] = true;
      //   }
      //   // the summation
      //   for (int i = 0; i < batch_size; i++) {
      //     if (tmp_flag == i+1) {
      //       pbbs::fetch_and_add(&overlaps[i], 1);
      //     }
      //   }
      // } // end parallel_for
      // size_t overlap_size = 0;
      // size_t overlap_size_one = 0;
      // size_t overlap_size_only = 0;
      // parallel_for(size_t j = 0; j < n; j++) {
      //   if (overlap_set[j]) {
      //     pbbs::fetch_and_add(&overlap_size, 1);
      //   }
      //   if (overlap_set_one[j]) {
      //     pbbs::fetch_and_add(&overlap_size_one, 1);
      //   }
      //   if (overlap_set_only[j]) {
      //     pbbs::fetch_and_add(&overlap_size_only, 1);
      //   }
      // } // end parallel_for
      // accumulated_overlap += overlap_size;
      // accumulated_overlap_one += overlap_size_one;
      // accumulated_overlap_only += overlap_size_only;
      // double total_overlap_score = 0.0;
      // for (int i = 0; i < batch_size; i++) {
      //   total_overlap_score += (i+1) * overlaps[i] * 1.0;
      // }
      // total_overlap_score = total_overlap_score / totalActivated / batch_size;
      // overlap_scores.push_back(total_overlap_score);

      // affinity_tracking.push_back(make_pair(Frontier.size(), 1.0 * overlap_size / Frontier.size()));
      // affinity_tracking_one.push_back(make_pair(Frontier.size(), 1.0 * overlap_size_one / Frontier.size()));
      // affinity_tracking_only.push_back(make_pair(Frontier.size(), 1.0 * overlap_size_only / Frontier.size()));
      // pbbs::delete_array(overlap_set, n);
      // pbbs::delete_array(overlap_set_one, n);
      // pbbs::delete_array(overlap_set_only, n);
    }

    // mode: no_dense, remove_duplicates (for batch size > 1)
    vertexSubset output = edgeMap(G, Frontier, BFSLV_F(Levels, CurrActiveArray, NextActiveArray, batch_size), -1, no_dense|remove_duplicates);

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
    bool isFEmpty = Frontier.isEmpty();
    Frontier.toDense();
    bool* new_d = Frontier.d;
    Frontier.d = nullptr;
    for(long i = 0; i < batch_size; i++) {
      if (defer_vec[i] == iteration) {
        new_d[vecQueries[i]] = true;
        NextActiveArray[(IdxType)vecQueries[i] * (IdxType)batch_size + (IdxType)i] = true;
      } 
      // else {
      //   if (isFEmpty && defer_vec[i] > iteration) {
      //     new_d[vecQueries[i]] = true;
      //     NextActiveArray[(IdxType)vecQueries[i] * (IdxType)batch_size + (IdxType)i] = true;
      //   }
      // }
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

  // profiling
  if (should_profile) {
    // double total_affinity = 0.0;
    // double total_overlap_score = 0.0;
    // cout << "Delay peak activation: " << peak_activation << endl;
    // cout << "Delay peak iteration: " << peak_iter << endl;

    // cout << "Delay Total iterations: " << iteration << endl;
    // cout << "Delay Total activations: " << totalActivated << endl;
    // double affinity_sum = 0.0;
    // for (int i = 0; i < affinity_tracking.size(); i++) {
    //   affinity_sum += affinity_tracking[i].first * affinity_tracking[i].second;
    // }
    // total_affinity = affinity_sum/totalActivated;
    // cout << "Delay AND Total affinity: " << total_affinity << endl;

    // double affinity_sum_one = 0.0;
    // double total_affinity_one = 0.0;
    // for (int i = 0; i < affinity_tracking_one.size(); i++) {
    //   affinity_sum_one += affinity_tracking_one[i].first * affinity_tracking_one[i].second;
    // }
    // total_affinity_one = affinity_sum_one/totalActivated;
    // cout << "Delay OR Total affinity: " << total_affinity_one << endl;

    // double affinity_sum_only = 0.0;
    // double total_affinity_only = 0.0;
    // for (int i = 0; i < affinity_tracking_only.size(); i++) {
    //   affinity_sum_only += affinity_tracking_only[i].first * affinity_tracking_only[i].second;
    // }
    // total_affinity_only = affinity_sum_only/totalActivated;
    // cout << "Delay [Only] affinity: " << total_affinity_only << endl;

    // // double total_overlap_score = 0.0;
    // for (int i = 0; i < batch_size; i++) {
    //   total_overlap_score += (i+1) * overlaps[i] * 1.0;
    // }
    // total_overlap_score = total_overlap_score / totalActivated / batch_size;
    // cout << "Delay Total Overlap score: " << total_overlap_score << endl;
    
    // cout << "Delay Iteration's overlap scores: " << endl;
    // for (int i = 0; i < overlap_scores.size(); i++) {
    //     cout << overlap_scores[i] << " ";
    // }
    // cout << endl;
    // int maxElementIndex = std::max_element(overlap_scores.begin(),overlap_scores.end()) - overlap_scores.begin();
    // double maxElement = *std::max_element(overlap_scores.begin(), overlap_scores.end());
    // cout << "Delay Max overlap score and iteration: " << maxElementIndex+1 << " " << maxElement << endl;
  }

#ifdef OUTPUT 
  for (int i = 0; i < batch_size; i++) {
    long start = vecQueries[i];
    char outFileName[300];
    sprintf(outFileName, "BFSLV_delay_output_src%ld.%ld.%ld.out", start, edge_count, batch_size);
    FILE *fp;
    fp = fopen(outFileName, "w");
    for (long j = 0; j < n; j++)
      fprintf(fp, "%ld %d\n", j, Levels[j * batch_size + i]);
    fclose(fp);
  }
#endif

  Frontier.del();
  pbbs::delete_array(CurrActiveArray, totalNumVertices);
  pbbs::delete_array(NextActiveArray, totalNumVertices);
  pbbs::delete_array(overlaps, batch_size);
  pbbs::delete_array(Levels, totalNumVertices);
  return make_pair(totalActivated, 0);
}

template <class vertex>
pair<size_t, size_t> Compute_Delay_Skipping(graph<vertex>& G, std::vector<long> vecQueries, commandLine P, std::vector<int> defer_vec, bool should_profile) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size;
  uintE* Levels = pbbs::new_array<uintE>(totalNumVertices);
  bool* frontier = pbbs::new_array<bool>(n);
  parallel_for(size_t i = 0; i < n; i++) {
    frontier[i] = false;
  }
  // for delaying initialization
  for(long i = 0; i < batch_size; i++) {
    if (defer_vec[i] == 0 || defer_vec[i] > 1000) {
      frontier[vecQueries[i]] = true;
    }
  }

  parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
    Levels[i] = (uintE)MAXLEVEL;
  }
  for(long i = 0; i < batch_size; i++) {
    Levels[(IdxType)batch_size * (IdxType)vecQueries[i] + (IdxType)i] = 0;
  }

  vertexSubset Frontier(n, frontier);

  // for profiling
  long iteration = 0;
  size_t totalActivated = 0;

  while(!Frontier.isEmpty()){
    iteration++;
    totalActivated += Frontier.size();
    // cout << "iteration: " << Frontier.size() << ": " << totalActivated << endl;
    // mode: no_dense, remove_duplicates (for batch size > 1)
    vertexSubset output = edgeMap(G, Frontier, BFSLV_SKIP_F(Levels, batch_size), -1, no_dense|remove_duplicates);

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
    sprintf(outFileName, "BFSLV_delay_skipping_output_src%ld.%ld.%ld.out", start, edge_count, batch_size);
    FILE *fp;
    fp = fopen(outFileName, "w");
    for (long j = 0; j < n; j++)
      fprintf(fp, "%ld %d\n", j, Levels[j * batch_size + i]);
    fclose(fp);
  }
#endif

  Frontier.del();
  pbbs::delete_array(Levels, totalNumVertices);
  return make_pair(totalActivated, 0);
}