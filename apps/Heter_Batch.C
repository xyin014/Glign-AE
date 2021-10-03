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
            if (CurrActiveArray[s_begin + j] & benchmark_sssp) {
                // Extract original values resides in the value array
                intE sValue = QueryValues[s_begin + j];
                intE dValue = QueryValues[d_begin + j];
                // Relax the numeric value of vertex d
                intE newDist = sValue + edgeLen;  // edgeLen = 1
                if (dValue > newDist) {
                    QueryValues[d_begin + j] = newDist;
                    NextActiveArray[d_begin + j] = benchmark_sssp;
                    ret = true;
                }
            } else if (CurrActiveArray[s_begin + j] & benchmark_bfs) {
                intE sValue = QueryValues[s_begin + j];
                intE dValue = QueryValues[d_begin + j];
                // Relax the numeric value of vertex d
                intE newDist = sValue + 1;  // edgeLen = 1

                if (dValue > newDist) {
                    QueryValues[d_begin + j] = newDist;
                    NextActiveArray[d_begin + j] = benchmark_bfs;
                    ret = true;
                }
            } else if (CurrActiveArray[s_begin + j] & benchmark_sswp) {
                // Extract original values resides in the value array
                intE sValue = QueryValues[s_begin + j];
                intE dValue = QueryValues[d_begin + j];
                // Relax the numeric value of vertex d
                intE newWidth = std::min(sValue, edgeLen);
                // Propagation of value and case of hub vertex.
                if (dValue < newWidth) {
                    QueryValues[d_begin + j] = newWidth;
                    NextActiveArray[d_begin + j] = benchmark_sswp;
                    ret = true;                    
                }
            }  else if (CurrActiveArray[s_begin + j] & benchmark_ssnp) {
                // Extract original values resides in the value array
                intE sValue = QueryValues[s_begin + j];
                intE dValue = QueryValues[d_begin + j];
                // Relax the numeric value of vertex d
                intE newWidth = std::max(sValue, edgeLen);
                // Propagation of value and case of hub vertex.
                if (dValue > newWidth) {
                    QueryValues[d_begin + j] = newWidth;
                    NextActiveArray[d_begin + j] = benchmark_sswp;
                    ret = true;                    
                }
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
            if (cur_active & benchmark_sssp) {
                // Extract original values resides in the value array
                intE sValue = QueryValues[s_begin + j];
                // Relax the numeric value of vertex d
                intE newDist = sValue + edgeLen;  // edgeLen = 1
                if (writeMin(&QueryValues[d_begin + j], newDist)
                        && CAS(&NextActiveArray[d_begin + j], _empty, benchmark_sssp)) {
                    ret = true;
                }
            } else if (cur_active & benchmark_bfs) {
                    // Extract original values resides in the value array
                    intE sValue = QueryValues[s_begin + j];
                    // Relax the numeric value of vertex d
                    intE newDist = sValue + 1;  // edgeLen = 1
                    if (writeMin(&QueryValues[d_begin + j], newDist)
                            && CAS(&NextActiveArray[d_begin + j], _empty, benchmark_bfs)) {
                        ret = true;
                    }
            } else if (cur_active & benchmark_sswp) {
                // Extract original values resides in the value array
                intE sValue = QueryValues[s_begin + j];
                // Relax the numeric value of vertex d
                intE newWidth = std::min(sValue, edgeLen);
                // Propagation of value and case of hub vertex.
                if (writeMax(&QueryValues[d_begin + j], newWidth)
                        && CAS(&NextActiveArray[d_begin + j], _empty, benchmark_sswp)) {
                    ret = true;
                }
            } else if (cur_active & benchmark_ssnp) {
                // Extract original values resides in the value array
                intE sValue = QueryValues[s_begin + j];
                // Relax the numeric value of vertex d
                intE newWidth = std::max(sValue, edgeLen);
                // Propagation of value and case of hub vertex.
                if (writeMin(&QueryValues[d_begin + j], newWidth)
                        && CAS(&NextActiveArray[d_begin + j], _empty, benchmark_sswp)) {
                    ret = true;
                }
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
uintE* Compute_Eval_Prop(graph<vertex>& G, std::vector<long> vecQueries, commandLine P) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size;

  uintE* ret = pbbs::new_array<uintE>(totalNumVertices);
  return ret;
}

template <class vertex>
void Compute_Base(graph<vertex>& G, std::vector<long> vecQueries, commandLine P, bool should_profile) {

}

template <class vertex>
void Compute_Heter_Base(graph<vertex>& G, std::vector<pair<long, benchmarkType>> vecQueries, commandLine P, bool should_profile) {
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
      if (Frontier.size() > peak_activation) {
        peak_activation = Frontier.size();
        peak_iter = iteration;
      }

      bool* overlap_set = pbbs::new_array<bool>(n); // activated for all queries
      bool* overlap_set_one = pbbs::new_array<bool>(n); // only activated for at least half of queries
      bool* overlap_set_only = pbbs::new_array<bool>(n);
      parallel_for(size_t i = 0; i < n; i++) {
        overlap_set[i] = true;
        overlap_set_one[i] = false;
        overlap_set_only[i] = false;
      }
      parallel_for(size_t index = 0; index < n; index++) {
        int tmp_flag = 0;
        for (int i = 0; i < batch_size; i++) {
          // if the vertex is activated for all queries...
          overlap_set[index] = overlap_set[index] && CurrActiveArray[index * batch_size + i];
          if (CurrActiveArray[index * batch_size + i]) {
            tmp_flag++;
          }
        }
        // activated for at least half of the queries.
        if (tmp_flag >= batch_size/2) {
          overlap_set_one[index] = true;
        }
        // only one is activated.
        if (tmp_flag == 1) {
          overlap_set_only[index] = true;
        }
        // the summation
        for (int i = 0; i < batch_size; i++) {
          if (tmp_flag == i+1) {
            pbbs::fetch_and_add(&overlaps[i], 1);
          }
        }
      } // end parallel_for
      size_t overlap_size = 0;
      size_t overlap_size_one = 0;
      size_t overlap_size_only = 0;
      parallel_for(size_t j = 0; j < n; j++) {
        if (overlap_set[j]) {
          pbbs::fetch_and_add(&overlap_size, 1);
        }
        if (overlap_set_one[j]) {
          pbbs::fetch_and_add(&overlap_size_one, 1);
        }
        if (overlap_set_only[j]) {
          pbbs::fetch_and_add(&overlap_size_only, 1);
        }
      } // end parallel_for
      accumulated_overlap += overlap_size;
      accumulated_overlap_one += overlap_size_one;
      accumulated_overlap_only += overlap_size_only;
      double total_overlap_score = 0.0;
      for (int i = 0; i < batch_size; i++) {
        total_overlap_score += (i+1) * overlaps[i] * 1.0;
      }
      total_overlap_score = total_overlap_score / totalActivated / batch_size;
      overlap_scores.push_back(total_overlap_score);

      affinity_tracking.push_back(make_pair(Frontier.size(), 1.0 * overlap_size / Frontier.size()));
      affinity_tracking_one.push_back(make_pair(Frontier.size(), 1.0 * overlap_size_one / Frontier.size()));
      affinity_tracking_only.push_back(make_pair(Frontier.size(), 1.0 * overlap_size_only / Frontier.size()));
      pbbs::delete_array(overlap_set, n);
      pbbs::delete_array(overlap_set_one, n);
      pbbs::delete_array(overlap_set_only, n);
    }

    // mode: no_dense, remove_duplicates (for batch size > 1)
    vertexSubset output = edgeMap(G, Frontier, Heterogeneous_F(QueryValues, CurrActiveArray, NextActiveArray, batch_size), -1, no_dense|remove_duplicates);

    Frontier.del();
    Frontier = output;

    // Frontier.toDense();
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
  if (should_profile) {
    double total_affinity = 0.0;
    double total_overlap_score = 0.0;
    cout << "Base peak activation: " << peak_activation << endl;
    cout << "Base peak iteration: " << peak_iter << endl;

    cout << "Base Total iterations: " << iteration << endl;
    cout << "Base Total activations: " << totalActivated << endl;
    double affinity_sum = 0.0;
    for (int i = 0; i < affinity_tracking.size(); i++) {
      affinity_sum += affinity_tracking[i].first * affinity_tracking[i].second;
    }
    total_affinity = affinity_sum/totalActivated;
    cout << "Base AND Total affinity: " << total_affinity << endl;

    double affinity_sum_one = 0.0;
    double total_affinity_one = 0.0;
    for (int i = 0; i < affinity_tracking_one.size(); i++) {
      affinity_sum_one += affinity_tracking_one[i].first * affinity_tracking_one[i].second;
    }
    total_affinity_one = affinity_sum_one/totalActivated;
    cout << "Base OR Total affinity: " << total_affinity_one << endl;

    double affinity_sum_only = 0.0;
    double total_affinity_only = 0.0;
    for (int i = 0; i < affinity_tracking_only.size(); i++) {
      affinity_sum_only += affinity_tracking_only[i].first * affinity_tracking_only[i].second;
    }
    total_affinity_only = affinity_sum_only/totalActivated;
    cout << "Base [Only] affinity: " << total_affinity_only << endl;

    // double total_overlap_score = 0.0;
    for (int i = 0; i < batch_size; i++) {
      total_overlap_score += (i+1) * overlaps[i] * 1.0;
    }
    total_overlap_score = total_overlap_score / totalActivated / batch_size;
    cout << "Base Total Overlap score: " << total_overlap_score << endl;
    
    cout << "Base Iteration's overlap scores: " << endl;
    for (int i = 0; i < overlap_scores.size(); i++) {
      cout << overlap_scores[i] << " ";
    }
    cout << endl;
    int maxElementIndex = std::max_element(overlap_scores.begin(),overlap_scores.end()) - overlap_scores.begin();
    double maxElement = *std::max_element(overlap_scores.begin(), overlap_scores.end());
    cout << "Base Max overlap score and iteration: " << maxElementIndex+1 << " " << maxElement << endl;

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
  pbbs::delete_array(QueryValues, totalNumVertices);
  pbbs::delete_array(CurrActiveArray, totalNumVertices);
  pbbs::delete_array(NextActiveArray, totalNumVertices);
  pbbs::delete_array(overlaps, batch_size);

}

template <class vertex>
void Compute_Delay(graph<vertex>& G, std::vector<long> vecQueries, commandLine P, std::vector<int> defer_vec, bool should_profile) {
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

    // profiling
    if (should_profile) {
      if (Frontier.size() > peak_activation) {
        peak_activation = Frontier.size();
        peak_iter = iteration;
      }
      bool* overlap_set = pbbs::new_array<bool>(n); // activated for all queries
      bool* overlap_set_one = pbbs::new_array<bool>(n); // only activated for at least half of queries
      bool* overlap_set_only = pbbs::new_array<bool>(n);
      parallel_for(size_t i = 0; i < n; i++) {
        overlap_set[i] = true;
        overlap_set_one[i] = false;
        overlap_set_only[i] = false;
      }
      parallel_for(size_t index = 0; index < n; index++) {
        int tmp_flag = 0;
        for (int i = 0; i < batch_size; i++) {
          // if the vertex is activated for all queries...
          overlap_set[index] = overlap_set[index] && CurrActiveArray[index * batch_size + i];
          if (CurrActiveArray[index * batch_size + i]) {
            tmp_flag++;
          }
        }
        // activated for at least half of the queries.
        if (tmp_flag >= batch_size/2) {
          overlap_set_one[index] = true;
        }
        // only one is activated.
        if (tmp_flag == 1) {
          overlap_set_only[index] = true;
        }
        // the summation
        for (int i = 0; i < batch_size; i++) {
          if (tmp_flag == i+1) {
            pbbs::fetch_and_add(&overlaps[i], 1);
          }
        }
      } // end parallel_for
      size_t overlap_size = 0;
      size_t overlap_size_one = 0;
      size_t overlap_size_only = 0;
      parallel_for(size_t j = 0; j < n; j++) {
        if (overlap_set[j]) {
          pbbs::fetch_and_add(&overlap_size, 1);
        }
        if (overlap_set_one[j]) {
          pbbs::fetch_and_add(&overlap_size_one, 1);
        }
        if (overlap_set_only[j]) {
          pbbs::fetch_and_add(&overlap_size_only, 1);
        }
      } // end parallel_for
      accumulated_overlap += overlap_size;
      accumulated_overlap_one += overlap_size_one;
      accumulated_overlap_only += overlap_size_only;
      double total_overlap_score = 0.0;
      for (int i = 0; i < batch_size; i++) {
        total_overlap_score += (i+1) * overlaps[i] * 1.0;
      }
      total_overlap_score = total_overlap_score / totalActivated / batch_size;
      overlap_scores.push_back(total_overlap_score);

      affinity_tracking.push_back(make_pair(Frontier.size(), 1.0 * overlap_size / Frontier.size()));
      affinity_tracking_one.push_back(make_pair(Frontier.size(), 1.0 * overlap_size_one / Frontier.size()));
      affinity_tracking_only.push_back(make_pair(Frontier.size(), 1.0 * overlap_size_only / Frontier.size()));
      pbbs::delete_array(overlap_set, n);
      pbbs::delete_array(overlap_set_one, n);
      pbbs::delete_array(overlap_set_only, n);

      // frontier_iterations.push_back(Frontier.size());
      // overlapped_iterations.push_back(overlap_size);
      // accumulated_overlapped_iterations.push_back(accumulated_overlap);
      // total_activated_iterations.push_back(totalActivated);
    }

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

  // profiling
  if (should_profile) {
    double total_affinity = 0.0;
    double total_overlap_score = 0.0;
    cout << "Delay peak activation: " << peak_activation << endl;
    cout << "Delay peak iteration: " << peak_iter << endl;

    cout << "Delay Total iterations: " << iteration << endl;
    cout << "Delay Total activations: " << totalActivated << endl;
    double affinity_sum = 0.0;
    for (int i = 0; i < affinity_tracking.size(); i++) {
      affinity_sum += affinity_tracking[i].first * affinity_tracking[i].second;
    }
    total_affinity = affinity_sum/totalActivated;
    cout << "Delay AND Total affinity: " << total_affinity << endl;

    double affinity_sum_one = 0.0;
    double total_affinity_one = 0.0;
    for (int i = 0; i < affinity_tracking_one.size(); i++) {
      affinity_sum_one += affinity_tracking_one[i].first * affinity_tracking_one[i].second;
    }
    total_affinity_one = affinity_sum_one/totalActivated;
    cout << "Delay OR Total affinity: " << total_affinity_one << endl;

    double affinity_sum_only = 0.0;
    double total_affinity_only = 0.0;
    for (int i = 0; i < affinity_tracking_only.size(); i++) {
      affinity_sum_only += affinity_tracking_only[i].first * affinity_tracking_only[i].second;
    }
    total_affinity_only = affinity_sum_only/totalActivated;
    cout << "Delay [Only] affinity: " << total_affinity_only << endl;

    // double total_overlap_score = 0.0;
    for (int i = 0; i < batch_size; i++) {
      total_overlap_score += (i+1) * overlaps[i] * 1.0;
    }
    total_overlap_score = total_overlap_score / totalActivated / batch_size;
    cout << "Delay Total Overlap score: " << total_overlap_score << endl;
    
    cout << "Delay Iteration's overlap scores: " << endl;
    for (int i = 0; i < overlap_scores.size(); i++) {
        cout << overlap_scores[i] << " ";
    }
    cout << endl;
    int maxElementIndex = std::max_element(overlap_scores.begin(),overlap_scores.end()) - overlap_scores.begin();
    double maxElement = *std::max_element(overlap_scores.begin(), overlap_scores.end());
    cout << "Delay Max overlap score and iteration: " << maxElementIndex+1 << " " << maxElement << endl;
  
    // // todo: remove
    // // vector<long> frontier_iterations;
    // // vector<long> overlapped_iterations;
    // // vector<long> accumulated_overlapped_iterations;
    // // vector<long> total_activated_iterations;
    // cout << "Delay frontier size per iteration: " << endl;
    // for (int i = 0; i < frontier_iterations.size(); i++) {
    //   cout << frontier_iterations[i] << " ";
    // }
    // cout << endl;

    // cout << "Delay overlapped size per iteration: " << endl;
    // for (int i = 0; i < overlapped_iterations.size(); i++) {
    //   cout << overlapped_iterations[i] << " ";
    // }
    // cout << endl;

    // cout << "Delay accumulated overlapped size per iteration: " << endl;
    // for (int i = 0; i < accumulated_overlapped_iterations.size(); i++) {
    //   cout << accumulated_overlapped_iterations[i] << " ";
    // }
    // cout << endl;

    // cout << "Delay total activated size per iteration: " << endl;
    // for (int i = 0; i < total_activated_iterations.size(); i++) {
    //   cout << total_activated_iterations[i] << " ";
    // }
    // cout << endl;
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
  pbbs::delete_array(overlaps, batch_size);
  pbbs::delete_array(ShortestPathLen, totalNumVertices);
}