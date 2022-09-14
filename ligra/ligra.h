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
#ifndef LIGRA_H
#define LIGRA_H
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>
#include <string>
#include <algorithm>
#include <vector>
#include <chrono>
#include <thread>
#include <cmath>
#include <random>
#include <queue>
#include <map>
#include <set>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>

#include "parallel.h"
#include "gettime.h"
#include "utils.h"
#include "vertex.h"
#include "compressedVertex.h"
#include "vertexSubset.h"
#include "graph.h"
#include "IO.h"
#include "parseCommandLine.h"
#include "index_map.h"
#include "edgeMap_utils.h"
using namespace std;

//*****START FRAMEWORK*****
#define MAXPATH 1000000
#define MAXLEVEL 1000000
#define MAXWIDTH 1000000

typedef uint32_t flags;
const flags no_output = 1;
const flags pack_edges = 2;
const flags sparse_no_filter = 4;
const flags dense_forward = 8;
const flags dense_parallel = 16;
const flags remove_duplicates = 32;
const flags no_dense = 64;
const flags edge_parallel = 128;
inline bool should_output(const flags& fl) { return !(fl & no_output); }

template <class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapDense(graph<vertex> GA, VS& vertexSubset, F &f, const flags fl) {
  using D = tuple<bool, data>;
  long n = GA.n;
  vertex *G = GA.V;
  if (should_output(fl)) {
    D* next = newA(D, n);
    auto g = get_emdense_gen<data>(next);
    parallel_for (long v=0; v<n; v++) {
      std::get<0>(next[v]) = 0;
      if (f.cond(v)) {
        G[v].decodeInNghBreakEarly(v, vertexSubset, f, g, fl & dense_parallel);
      }
    }
    return vertexSubsetData<data>(n, next);
  } else {
    auto g = get_emdense_nooutput_gen<data>();
    parallel_for (long v=0; v<n; v++) {
      if (f.cond(v)) {
        G[v].decodeInNghBreakEarly(v, vertexSubset, f, g, fl & dense_parallel);
      }
    }
    return vertexSubsetData<data>(n);
  }
}

template <class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapDenseForward(graph<vertex> GA, VS& vertexSubset, F &f, const flags fl) {
  using D = tuple<bool, data>;
  long n = GA.n;
  vertex *G = GA.V;
  if (should_output(fl)) {
    D* next = newA(D, n);
    auto g = get_emdense_forward_gen<data>(next);
    parallel_for(long i=0;i<n;i++) { std::get<0>(next[i]) = 0; }
    parallel_for (long i=0; i<n; i++) {
      if (vertexSubset.isIn(i)) {
        G[i].decodeOutNgh(i, f, g);
      }
    }
    return vertexSubsetData<data>(n, next);
  } else {
    auto g = get_emdense_forward_nooutput_gen<data>();
    parallel_for (long i=0; i<n; i++) {
      if (vertexSubset.isIn(i)) {
        G[i].decodeOutNgh(i, f, g);
      }
    }
    return vertexSubsetData<data>(n);
  }
}

template <class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapSparse(graph<vertex>& GA, vertex* frontierVertices, VS& indices,
        uintT* degrees, uintT m, F &f, const flags fl) {
  using S = tuple<uintE, data>;
  long n = indices.n;
  S* outEdges;
  long outEdgeCount = 0;

  if (should_output(fl)) {
    uintT* offsets = degrees;
    outEdgeCount = sequence::plusScan(offsets, offsets, m);
    outEdges = newA(S, outEdgeCount);
    auto g = get_emsparse_gen<data>(outEdges);
    parallel_for (size_t i = 0; i < m; i++) {
      uintT v = indices.vtx(i), o = offsets[i];
      vertex vert = frontierVertices[i];
      vert.decodeOutNghSparse(v, o, f, g);
    }
  } else {
    auto g = get_emsparse_nooutput_gen<data>();
    parallel_for (size_t i = 0; i < m; i++) {
      uintT v = indices.vtx(i);
      vertex vert = frontierVertices[i];
      vert.decodeOutNghSparse(v, 0, f, g);
    }
  }

  if (should_output(fl)) {
    S* nextIndices = newA(S, outEdgeCount);
    if (fl & remove_duplicates) {
      if (GA.flags == NULL) {
        GA.flags = newA(uintE, n);
        parallel_for(long i=0;i<n;i++) { GA.flags[i]=UINT_E_MAX; }
      }
      auto get_key = [&] (size_t i) -> uintE& { return std::get<0>(outEdges[i]); };
      remDuplicates(get_key, GA.flags, outEdgeCount, n);
    }
    auto p = [] (tuple<uintE, data>& v) { return std::get<0>(v) != UINT_E_MAX; };
    size_t nextM = pbbs::filterf(outEdges, nextIndices, outEdgeCount, p);
    free(outEdges);
    return vertexSubsetData<data>(n, nextM, nextIndices);
  } else {
    return vertexSubsetData<data>(n);
  }
}

template <class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapSparse_no_filter(graph<vertex>& GA,
    vertex* frontierVertices, VS& indices, uintT* offsets, uintT m, F& f,
    const flags fl) {
  using S = tuple<uintE, data>;
  long n = indices.n;
  long outEdgeCount = sequence::plusScan(offsets, offsets, m);
  S* outEdges = newA(S, outEdgeCount);

  auto g = get_emsparse_no_filter_gen<data>(outEdges);

  // binary-search into scan to map workers->chunks
  size_t b_size = 10000;
  size_t n_blocks = nblocks(outEdgeCount, b_size);

  uintE* cts = newA(uintE, n_blocks+1);
  size_t* block_offs = newA(size_t, n_blocks+1);

  auto offsets_m = make_in_imap<uintT>(m, [&] (size_t i) { return offsets[i]; });
  auto lt = [] (const uintT& l, const uintT& r) { return l < r; };
  parallel_for(size_t i=0; i<n_blocks; i++) {
    size_t s_val = i*b_size;
    block_offs[i] = pbbs::binary_search(offsets_m, s_val, lt);
  }
  block_offs[n_blocks] = m;
  parallel_for (size_t i=0; i<n_blocks; i++) {
    if ((i == n_blocks-1) || block_offs[i] != block_offs[i+1]) {
      // start and end are offsets in [m]
      size_t start = block_offs[i];
      size_t end = block_offs[i+1];
      uintT start_o = offsets[start];
      uintT k = start_o;
      for (size_t j=start; j<end; j++) {
        uintE v = indices.vtx(j);
        size_t num_in = frontierVertices[j].decodeOutNghSparseSeq(v, k, f, g);
        k += num_in;
      }
      cts[i] = (k - start_o);
    } else {
      cts[i] = 0;
    }
  }

  long outSize = sequence::plusScan(cts, cts, n_blocks);
  cts[n_blocks] = outSize;

  S* out = newA(S, outSize);

  parallel_for (size_t i=0; i<n_blocks; i++) {
    if ((i == n_blocks-1) || block_offs[i] != block_offs[i+1]) {
      size_t start = block_offs[i];
      size_t start_o = offsets[start];
      size_t out_off = cts[i];
      size_t block_size = cts[i+1] - out_off;
      for (size_t j=0; j<block_size; j++) {
        out[out_off + j] = outEdges[start_o + j];
      }
    }
  }
  free(outEdges); free(cts); free(block_offs);

  if (fl & remove_duplicates) {
    if (GA.flags == NULL) {
      GA.flags = newA(uintE, n);
      parallel_for(size_t i=0;i<n;i++) { GA.flags[i]=UINT_E_MAX; }
    }
    auto get_key = [&] (size_t i) -> uintE& { return std::get<0>(out[i]); };
    remDuplicates(get_key, GA.flags, outSize, n);
    S* nextIndices = newA(S, outSize);
    auto p = [] (tuple<uintE, data>& v) { return std::get<0>(v) != UINT_E_MAX; };
    size_t nextM = pbbs::filterf(out, nextIndices, outSize, p);
    free(out);
    return vertexSubsetData<data>(n, nextM, nextIndices);
  }
  return vertexSubsetData<data>(n, outSize, out);
}

// Decides on sparse or dense base on number of nonzeros in the active vertices.
template <class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapData(graph<vertex>& GA, VS &vs, F f,
    intT threshold = -1, const flags& fl=0) {
  long numVertices = GA.n, numEdges = GA.m, m = vs.numNonzeros();
  if(threshold == -1) threshold = numEdges/20; //default threshold
  vertex *G = GA.V;
  if (numVertices != vs.numRows()) {
    cout << "edgeMap: Sizes Don't match" << endl;
    abort();
  }
  if (m == 0) return vertexSubsetData<data>(numVertices);
  uintT* degrees = NULL;
  vertex* frontierVertices = NULL;
  uintT outDegrees = 0;
  if(threshold > 0) { //compute sum of out-degrees if threshold > 0 
    vs.toSparse();
    degrees = newA(uintT, m);
    frontierVertices = newA(vertex,m);
    {parallel_for (size_t i=0; i < m; i++) {
	uintE v_id = vs.vtx(i);
	vertex v = G[v_id];
	degrees[i] = v.getOutDegree();
	frontierVertices[i] = v;
      }}
    outDegrees = sequence::plusReduce(degrees, m);
    if (outDegrees == 0) return vertexSubsetData<data>(numVertices);
  }
   if (!(fl & no_dense) && m + outDegrees > threshold) {
    // cout << "dense mod\n";
    if(degrees) free(degrees);
    if(frontierVertices) free(frontierVertices);
    vs.toDense();
    return (fl & dense_forward) ?
      edgeMapDenseForward<data, vertex, VS, F>(GA, vs, f, fl) :
      edgeMapDense<data, vertex, VS, F>(GA, vs, f, fl);
   } 
   else {
     auto vs_out =
       (should_output(fl) && fl & sparse_no_filter) ? // only call snof when we output
       edgeMapSparse_no_filter<data, vertex, VS, F>(GA, frontierVertices, vs, degrees, vs.numNonzeros(), f, fl) :
       edgeMapSparse<data, vertex, VS, F>(GA, frontierVertices, vs, degrees, vs.numNonzeros(), f, fl);
     free(degrees); free(frontierVertices);
     return vs_out;
   }
}

// Regular edgeMap, where no extra data is stored per vertex.
template <class vertex, class VS, class F>
vertexSubset edgeMap(graph<vertex>& GA, VS& vs, F f,
    intT threshold = -1, const flags& fl=0) {
  return edgeMapData<pbbs::empty>(GA, vs, f, threshold, fl);
}

// Packs out the adjacency lists of all vertex in vs. A neighbor, ngh, is kept
// in the new adjacency list if p(ngh) is true.
// Weighted graphs are not yet supported, but this should be easy to do.
template <class vertex, class P>
vertexSubsetData<uintE> packEdges(graph<vertex>& GA, vertexSubset& vs, P& p, const flags& fl=0) {
  using S = tuple<uintE, uintE>;
  vs.toSparse();
  vertex* G = GA.V; long m = vs.numNonzeros(); long n = vs.numRows();
  if (vs.size() == 0) {
    return vertexSubsetData<uintE>(n);
  }
  auto degrees = array_imap<uintT>(m);
  granular_for(i, 0, m, (m > 2000), {
    uintE v = vs.vtx(i);
    degrees[i] = G[v].getOutDegree();
  });
  long outEdgeCount = pbbs::scan_add(degrees, degrees);
  S* outV;
  if (should_output(fl)) {
    outV = newA(S, vs.size());
  }

  bool* bits = newA(bool, outEdgeCount);
  uintE* tmp1 = newA(uintE, outEdgeCount);
  uintE* tmp2 = newA(uintE, outEdgeCount);
  if (should_output(fl)) {
    parallel_for (size_t i=0; i<m; i++) {
      uintE v = vs.vtx(i);
      size_t offset = degrees[i];
      auto bitsOff = &(bits[offset]); auto tmp1Off = &(tmp1[offset]);
      auto tmp2Off = &(tmp2[offset]);
      size_t ct = G[v].packOutNgh(v, p, bitsOff, tmp1Off, tmp2Off);
      outV[i] = make_tuple(v, ct);
    }
  } else {
    parallel_for (size_t i=0; i<m; i++) {
      uintE v = vs.vtx(i);
      size_t offset = degrees[i];
      auto bitsOff = &(bits[offset]); auto tmp1Off = &(tmp1[offset]);
      auto tmp2Off = &(tmp2[offset]);
      size_t ct = G[v].packOutNgh(v, p, bitsOff, tmp1Off, tmp2Off);
    }
  }
  free(bits); free(tmp1); free(tmp2);
  if (should_output(fl)) {
    return vertexSubsetData<uintE>(n, m, outV);
  } else {
    return vertexSubsetData<uintE>(n);
  }
}

template <class vertex, class P>
vertexSubsetData<uintE> edgeMapFilter(graph<vertex>& GA, vertexSubset& vs, P& p, const flags& fl=0) {
  vs.toSparse();
  if (fl & pack_edges) {
    return packEdges<vertex, P>(GA, vs, p, fl);
  }
  vertex* G = GA.V; long m = vs.numNonzeros(); long n = vs.numRows();
  using S = tuple<uintE, uintE>;
  if (vs.size() == 0) {
    return vertexSubsetData<uintE>(n);
  }
  S* outV;
  if (should_output(fl)) {
    outV = newA(S, vs.size());
  }
  if (should_output(fl)) {
    parallel_for (size_t i=0; i<m; i++) {
      uintE v = vs.vtx(i);
      size_t ct = G[v].countOutNgh(v, p);
      outV[i] = make_tuple(v, ct);
    }
  } else {
    parallel_for (size_t i=0; i<m; i++) {
      uintE v = vs.vtx(i);
      size_t ct = G[v].countOutNgh(v, p);
    }
  }
  if (should_output(fl)) {
    return vertexSubsetData<uintE>(n, m, outV);
  } else {
    return vertexSubsetData<uintE>(n);
  }
}



//*****VERTEX FUNCTIONS*****

template <class F, class VS, typename std::enable_if<
  !std::is_same<VS, vertexSubset>::value, int>::type=0 >
void vertexMap(VS& V, F f) {
  size_t n = V.numRows(), m = V.numNonzeros();
  if(V.dense()) {
    parallel_for(long i=0;i<n;i++) {
      if(V.isIn(i)) {
        f(i, V.ithData(i));
      }
    }
  } else {
    parallel_for(long i=0;i<m;i++) {
      f(V.vtx(i), V.vtxData(i));
    }
  }
}

template <class VS, class F, typename std::enable_if<
  std::is_same<VS, vertexSubset>::value, int>::type=0 >
void vertexMap(VS& V, F f) {
  size_t n = V.numRows(), m = V.numNonzeros();
  if(V.dense()) {
    parallel_for(long i=0;i<n;i++) {
      if(V.isIn(i)) {
        f(i);
      }
    }
  } else {
    parallel_for(long i=0;i<m;i++) {
      f(V.vtx(i));
    }
  }
}

//Note: this is the version of vertexMap in which only a subset of the
//input vertexSubset is returned
template <class F>
vertexSubset vertexFilter(vertexSubset V, F filter) {
  long n = V.numRows(), m = V.numNonzeros();
  V.toDense();
  bool* d_out = newA(bool,n);
  {parallel_for(long i=0;i<n;i++) d_out[i] = 0;}
  {parallel_for(long i=0;i<n;i++)
      if(V.d[i]) d_out[i] = filter(i);}
  return vertexSubset(n,d_out);
}

template <class F>
vertexSubset vertexFilter2(vertexSubset V, F filter) {
  long n = V.numRows(), m = V.numNonzeros();
  if (m == 0) {
    return vertexSubset(n);
  }
  bool* bits = newA(bool, m);
  V.toSparse();
  {parallel_for(size_t i=0; i<m; i++) {
    uintE v = V.vtx(i);
    bits[i] = filter(v);
  }}
  auto v_imap = make_in_imap<uintE>(m, [&] (size_t i) { return V.vtx(i); });
  auto bits_m = make_in_imap<bool>(m, [&] (size_t i) { return bits[i]; });
  auto out = pbbs::pack(v_imap, bits_m);
  out.alloc = false;
  free(bits);
  return vertexSubset(n, out.size(), out.s);
}

template <class data, class F>
vertexSubset vertexFilter2(vertexSubsetData<data> V, F filter) {
  long n = V.numRows(), m = V.numNonzeros();
  if (m == 0) {
    return vertexSubset(n);
  }
  bool* bits = newA(bool, m);
  V.toSparse();
  parallel_for(size_t i=0; i<m; i++) {
    auto t = V.vtxAndData(i);
    bits[i] = filter(std::get<0>(t), std::get<1>(t));
  }
  auto v_imap = make_in_imap<uintE>(m, [&] (size_t i) { return V.vtx(i); });
  auto bits_m = make_in_imap<bool>(m, [&] (size_t i) { return bits[i]; });
  auto out = pbbs::pack(v_imap, bits_m);
  out.alloc = false;
  free(bits);
  return vertexSubset(n, out.size(), out.s);
}



//cond function that always returns true
inline bool cond_true (intT d) { return 1; }

template<class vertex>
void Compute(graph<vertex>&, commandLine);

template<class vertex>
uintE* Compute_Eval(graph<vertex>&, vector<long>, commandLine);
// template<class vertex>
// uintE* Compute_Eval_Prop(graph<vertex>&, vector<long>, commandLine);
template<class vertex>
pair<size_t, size_t> Compute_Base(graph<vertex>&, vector<long>, commandLine, bool should_profile=false);
template<class vertex>
pair<size_t, size_t> Compute_Chunk(graph<vertex>&, vector<long>, set<set<long>>&, commandLine, bool should_profile=false);
template<class vertex>
pair<size_t, size_t> Compute_Chunk_v2(graph<vertex>&, vector<long>, set<set<long>>&, commandLine, bool should_profile=false);
template<class vertex>
pair<size_t, size_t> Compute_Chunk_v3(graph<vertex>&, vector<long>, set<set<long>>&, commandLine, bool should_profile=false);
template<class vertex>
pair<size_t, size_t> Compute_Delay(graph<vertex>&, vector<long>, commandLine, vector<int>, bool should_profile=false);
template<class vertex>
pair<size_t, size_t> Compute_Delay_Skipping(graph<vertex>&, vector<long>, commandLine, vector<int>, bool should_profile=false);
template<class vertex>
pair<size_t, size_t> Compute_Base_Skipping(graph<vertex>&, vector<long>, commandLine, bool should_profile=false);
template<class vertex>
pair<size_t, size_t> Compute_Base_Dynamic(graph<vertex>&, vector<long>, queue<long>&, commandLine, bool should_profile=false);
template<class vertex>
pair<size_t, size_t> Compute_Chunk_Concurrent(graph<vertex>&, vector<long>&, intE*, bool*, bool*, long, long, commandLine, bool should_profile=false);

template<class vertex>
pair<vector<long>, vector<long>> streamingPreprocessingReturnHops(graph<vertex>& G, vector<long> userQueries, int n_high_deg, int combination_max, commandLine P, uintE* distances);

template<class vertex>
void Compute(hypergraph<vertex>&, commandLine);

bool sortByLargerSecondElement(const pair<long, long> &a, const pair<long, long> &b) {
  return (a.second > b.second);
}

void groundtruth(int argc, char* argv[]) {
  commandLine P(argc,argv," [-s] <inFile>");
  char* iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");
  bool compressed = P.getOptionValue("-c");
  bool binary = P.getOptionValue("-b");
  bool mmap = P.getOptionValue("-m");
  //cout << "mmap = " << mmap << endl;
  long rounds = P.getOptionLongValue("-rounds",1);
  long max_delay = P.getOptionLongValue("-max_delay", 8);

  string queryFileName = string(P.getOptionValue("-qf", ""));
  int combination_max = P.getOptionLongValue("-max_combination", 256);
  size_t bSize = P.getOptionLongValue("-batch", 4);
  size_t n_high_deg = P.getOptionLongValue("-nhighdeg", 4);

  cout << "graph file name: " << iFile << endl;
  cout << "query file name: " << queryFileName << endl;

  // Initialization and preprocessing
  std::vector<long> userQueries; 
  long start = -1;
  char inFileName[300];
  ifstream inFile;
  sprintf(inFileName, queryFileName.c_str());
  inFile.open(inFileName, ios::in);
  while (inFile >> start) {
    userQueries.push_back(start);
  }
  inFile.close();

  std::random_device rd;
  auto rng = std::default_random_engine { rd() };

  cout << "number of random queries: " << userQueries.size() << endl;
  int batch_size = userQueries.size();
  
  if (symmetric) {
    cout << "symmetric graph: not implemented!\n";

  } else {
    // For directed graph...
    cout << "asymmetric graph\n";
    graph<asymmetricVertex> G =
      readGraph<asymmetricVertex>(iFile,compressed,symmetric,binary,mmap); //asymmetric graph
    cout << "n=" << G.n << " m=" << G.m << endl;

    size_t n = G.n;
    // ========================================
    // finding the high degree vertices and evluating BFS on high degree vtxs
    std::vector<std::pair<long, long>> vIDDegreePairs;
    for (long i = 0; i < n; i++) {
      long temp_degree =  G.V[i].getOutDegree();
      if (temp_degree >= 50) {
        vIDDegreePairs.push_back(std::make_pair(i, temp_degree));
      }
    }
    std::sort(vIDDegreePairs.begin(), vIDDegreePairs.end(), sortByLargerSecondElement);
    vector<long> highdegQ;
    int high_deg_batch = n_high_deg;
    cout << "High Deg Vtxs: \n";
    for (int i = 0; i < high_deg_batch; i++) {
      highdegQ.push_back(vIDDegreePairs[i].first);
      cout << vIDDegreePairs[i].first << ", degree: " << vIDDegreePairs[i].second << endl;
    }

    uintE* distances_multiple;
    // On edge reversed graph...
    G.transpose();
    distances_multiple = Compute_Eval(G,highdegQ,P);  // to get hops
    uintE* distances = pbbs::new_array<uintE>(n);
    G.transpose();
    parallel_for(size_t i = 0; i < n; i++) {
      distances[i] = (uintE)MAXLEVEL;
    }
    parallel_for(size_t i = 0; i < n; i++) {
      for (int j = 0; j < high_deg_batch; j++) {
        if (distances_multiple[j+i*high_deg_batch] < distances[i]) {
          distances[i] = distances_multiple[j+i*high_deg_batch];
        }
      }
    }
    // hop distributions of input queries.
    std::map<long, long> user_hist;
    for (long i = 0; i < userQueries.size(); i++) {
      int dist = distances[userQueries[i]];
      user_hist[dist]++;
    }
    for (const auto& x : user_hist) std::cout << x.first << " " << x.second <<"\n";
    
    // Query evaluation: sequential, batching, delayed batching
    vector<long> batchedQuery;
    batchedQuery = userQueries;
    cout << "=================\n";
    vector<double> diff_map_intra(16, 0); // first diff in Table 13.
    vector<double> diff_map_inter(16, 0);
    vector<double> diff_map_best(16, 0);
    map<int, int> diff_cnt;

    for (int i = 0; i < combination_max; i++) {
      timer t_seq, t_batch;
      std::shuffle(std::begin(batchedQuery), std::end(batchedQuery), rng);
      vector<long> tmp_batch;
      cout << "Evaluating queries: ";
      for (int j = 0; j < bSize; j++) {
        long tmp_query_id = batchedQuery[j];
        tmp_batch.push_back(tmp_query_id);
        cout << tmp_query_id << " ";
      }
      cout << endl;

      // Sequential
      long seq_cnt = 0;
      t_seq.start();
      for (int j = 0; j < tmp_batch.size(); j++) {
        vector<long> tmp_single_query;
        tmp_single_query.push_back(tmp_batch[j]);
        // Compute_Base(G,tmp_single_query,P);
        auto tmp_cnt = Compute_Base_Skipping(G,tmp_single_query,P,0);
        seq_cnt += tmp_cnt.first;
      }
      t_seq.stop();
      double seq_time = t_seq.totalTime;
      cout << "Sequential total activations: " << seq_cnt << endl;

      // Batching
      long batch_cnt = 0;
      t_batch.start();
      for (int r = 0; r < rounds; r++) {
        // Compute_Base(G,tmp_batch,P);
        auto tmp_cnt = Compute_Base_Skipping(G,tmp_batch,P,0);
        batch_cnt = tmp_cnt.first;
      }
      t_batch.stop();
      cout << "Baseline batching total activations: " << batch_cnt << endl;
      cout << "Baseline batching affinity: " << 1.0*batch_cnt / seq_cnt << endl;

      vector<int> tmp_delay = {0, 0};
      int max_ind = distances[tmp_batch[0]] > distances[tmp_batch[1]] ? 0 : 1;
      int min_ind = distances[tmp_batch[0]] < distances[tmp_batch[1]] ? 0 : 1;
      vector<double> combined_spd;
      vector<double> combined_affinity;
      for (int j = 0; j < max_delay; j++) {
        tmp_delay[min_ind] = j;
        cout << "No. " << min_ind << " delays " << j << " iterations\n";
        long delay_cnt = 0;
        timer t_delay; t_delay.start();
        // Compute_Delay(G,tmp_batch,P,tmp_delay);
        auto tmp_cnt = Compute_Delay_Skipping(G,tmp_batch,P,tmp_delay);

        // Compute_Base(G,tmp_batch,P);
        t_delay.stop();
        delay_cnt = tmp_cnt.first;
        double delay_time = t_delay.totalTime;
        combined_spd.push_back(seq_time / delay_time);
        combined_affinity.push_back(1.0*delay_cnt / seq_cnt);
        cout << "delayed query time: " << delay_time << endl;
        cout << "Delayed Eval Speedup: " << seq_time / delay_time << endl;
      }

      // // Delayed batching
      vector<int> dist_to_high;
      long total_delays = 0;
      for (int j = 0; j < tmp_batch.size(); j++) {
        cout << "q" << j << " to highest deg vtx: " << distances[tmp_batch[j]] << endl;
        if (distances[tmp_batch[j]] != MAXLEVEL) {
          dist_to_high.push_back(distances[tmp_batch[j]]);
        } else {
          dist_to_high.push_back(-1);
        }
      }
      int max_dist_to_high = *max_element(dist_to_high.begin(), dist_to_high.end());

      for (int j = 0; j < dist_to_high.size(); j++) {
        if (dist_to_high[j] == -1) {
          dist_to_high[j] = max_dist_to_high;
        }
      }

      for (int j = 0; j < dist_to_high.size(); j++) {
        dist_to_high[j] = max_dist_to_high - dist_to_high[j];
        total_delays += dist_to_high[j];
      }
      
      double batch_time = t_batch.totalTime / rounds;
      t_seq.reportTotal("sequential time");
      t_batch.reportTotal("batching evaluation time");

      cout << "Batching speedup: " << seq_time / batch_time << endl;

      // Print results
      cout << "Combined results for delaying: ";
      for (int j = 0; j < combined_spd.size(); j++) {
        cout << combined_spd[j] << " ";
      }
      cout << endl;
      cout << "Combined results for affinity: ";
      for (int j = 0; j < combined_spd.size(); j++) {
        cout << combined_affinity[j] << " ";
      }
      cout << endl;

      int best_delay = -1;
      double best_spd = 0.0;
      double best_affinity = 0.0;
      for (int j = 0; j < combined_spd.size(); j++) {
        if (combined_spd[j] > best_spd) {
          best_spd = combined_spd[j];
          best_affinity = combined_affinity[j];
          best_delay = j;
        }
      }
      int tmp_delta = abs(dist_to_high[0] - dist_to_high[1]);
      cout << "Best delay: " << best_delay << " " << best_spd  << " " << best_affinity << endl;
      cout << "Delta_Speedup: " << tmp_delta << " " << combined_spd[tmp_delta] << " " << combined_affinity[tmp_delta] << endl;
      cout << "=================\n";
      int diff = abs(best_delay-tmp_delta);
      diff_map_intra[diff] += combined_spd[0];
      diff_map_inter[diff] += combined_spd[tmp_delta];
      diff_map_best[diff] += best_spd;
      diff_cnt[diff]++;
    }
    
    cout << "Diff\t Cnt\t Ratio\t Glign-Intra\t Glign-Inter\t Best-Align\n";
    for (int i = 0; i < diff_map_intra.size(); i++) {
      cout << i << "\t " << diff_cnt[i] << "\t " << diff_map_intra[i]/diff_cnt[i] << "\t " << diff_map_inter[i]/diff_cnt[i] << "\t " << diff_map_best[i]/diff_cnt[i] << endl;
    }
    
    G.del();
    pbbs::delete_array(distances, n);
    pbbs::delete_array(distances_multiple, n*highdegQ.size());
  }
}

template <class vertex>
pair<vector<long>, vector<long>> streamingPreprocessing(graph<vertex>& G, vector<long> userQueries, int n_high_deg, int combination_max, commandLine P) {
  // input: graph, queries, queries_to_process, number of high degree vtxs.
  // return: unsorted and sorted queries of size queries_to_process
  vector<long> truncatedQueries;
  vector<long> sortedQueries;
  size_t n = G.n;
  std::vector<std::pair<long, long>> vIDDegreePairs;
  
  for (long i = 0; i < combination_max; i++) {
    truncatedQueries.push_back(userQueries[i]);
  }
  return make_pair(truncatedQueries, sortedQueries);
}

template <class vertex>
pair<vector<long>, vector<long>> streamingPreprocessingReturnHops(graph<vertex>& G, vector<long> userQueries, int n_high_deg, int combination_max, commandLine P, uintE* distances) {
  // input: graph, queries, queries_to_process, number of high degree vtxs.
  // return: unsorted and sorted queries of size |queries_to_process|
  vector<long> truncatedQueries;
  vector<long> sortedQueries;
  size_t n = G.n;
  std::vector<std::pair<long, long>> vIDDegreePairs;
  for (long i = 0; i < n; i++) {
    long temp_degree =  G.V[i].getOutDegree();
    if (temp_degree >= 1) {
      vIDDegreePairs.push_back(std::make_pair(i, temp_degree));
    }
  }
  std::sort(vIDDegreePairs.begin(), vIDDegreePairs.end(), sortByLargerSecondElement);
  vector<long> highdegQ;
  int high_deg_batch = n_high_deg;
  for (int i = 0; i < high_deg_batch; i++) {
    highdegQ.push_back(vIDDegreePairs[i].first);
  }

  uintE* distances_multiple;
  // On edge reversed graph...
  timer t_t1;
  t_t1.start();
  G.transpose();
  distances_multiple = Compute_Eval(G,highdegQ,P);  // to get hops (bfs distance)
  G.transpose();
  
  parallel_for(size_t i = 0; i < n; i++) {
    distances[i] = (uintE)MAXLEVEL;
  }
  parallel_for(size_t i = 0; i < n; i++) {
    for (int j = 0; j < high_deg_batch; j++) {
      if (distances_multiple[j+i*high_deg_batch] < distances[i]) {
        distances[i] = distances_multiple[j+i*high_deg_batch];
      }
    }
  }
  // hop distributions of input queries.
  std::map<long, long> user_hist;
  for (long i = 0; i < userQueries.size(); i++) {
    int dist = distances[userQueries[i]];
    user_hist[dist]++;
  }
  for (long i = 0; i < combination_max; i++) {
    truncatedQueries.push_back(userQueries[i]);
  }
  std::vector<std::pair<size_t, long>> vtxHopPairs;
  for (long i = 0; i < truncatedQueries.size(); i++) {
    vtxHopPairs.push_back(std::make_pair(truncatedQueries[i], distances[truncatedQueries[i]]));
  }
  std::sort(vtxHopPairs.begin(), vtxHopPairs.end(), sortByLargerSecondElement);
  for (int i = 0; i < vtxHopPairs.size(); i++) {
    long vID = vtxHopPairs[i].first;
    sortedQueries.push_back(vID);
  }
  t_t1.stop();
  double time1 = t_t1.totalTime;
  cout << "Profiling cost: " << time1 << endl;
  return make_pair(truncatedQueries, sortedQueries);
}

template <class vertex>
vector<pair<size_t, size_t>> bufferStreaming(graph<vertex>& G, std::vector<long> bufferedQueries, int bSize, commandLine P) {
  vector<pair<size_t, size_t>> res;
  timer t_t1;
  t_t1.start();
  for (int i = 0; i < bufferedQueries.size(); i=i+bSize) {
    std::vector<long> tmpBatch;
    cout << "evalating queries: " << endl;
    for (int j = 0; j < bSize; j++) {
      cout << bufferedQueries[i+j] << " ";
      tmpBatch.push_back(bufferedQueries[i+j]);
    }
    cout << endl;
    pair<size_t, size_t> share_cnt = Compute_Base(G,tmpBatch,P,true);
    res.push_back(share_cnt);
  }
  t_t1.stop();
  double time1 = t_t1.totalTime;
  cout << "Ligra-C time: " << time1 << endl;
  return res;
}

template <class vertex>
vector<pair<size_t, size_t>> streamingWithCGQ(graph<vertex>& G, std::vector<long> bufferedQueries, set<set<long>>& C_Set, vector<long>& chunk_lookup, int bSize, commandLine P, bool should_profile=false) {
  timer start_time1; start_time1.start();
  for (int i = 0; i < bufferedQueries.size(); i=i+bSize) {
    // cout << "i: " << i << endl;
    std::vector<long> tmpBatch;
    for (int j = 0; j < bSize; j++) {
      tmpBatch.push_back(bufferedQueries[i+j]);
    }
    
    timer t_t1;
    t_t1.start();
    if (bSize != 1)
      // Compute_Base(G,tmpBatch,P);
    t_t1.stop();
    double time1 = t_t1.totalTime;   
  }
  start_time1.stop();
  double query_time1 = start_time1.totalTime;

  vector<pair<size_t, size_t>> res;
  long n = G.n;
  timer t_t1;
  t_t1.start();
  std::vector<long> tmpBatch;
  for (int i = 0; i < bufferedQueries.size(); i=i+bSize) {
    for (int j = 0; j < bSize; j++) {
      tmpBatch.push_back(bufferedQueries[i+j]);
    }
    {parallel_for (int j = 0; j < bSize; j++) {
      std::vector<long> tmpBatch;
      tmpBatch.push_back(bufferedQueries[i+j]);
      Compute_Base_Skipping(G,tmpBatch, P,0);
    }}
    cout << "finished concurrent queries\n ";
  }
  t_t1.stop();
  double time1 = t_t1.totalTime;
  cout << "Profiling time: " << time1 << endl;
  return res;
}

template <class vertex>
vector<pair<size_t, size_t>> bufferStreamingSkipping(graph<vertex>& G, std::vector<long> bufferedQueries, int bSize, commandLine P, uintE* distances, bool should_profile=false) {
  bool shouldDelay = P.getOptionValue("-delay");
  long selection = P.getOptionLongValue("-mode",1);
  vector<pair<size_t, size_t>> res;
  double batching_time = 0;
  for (int i = 0; i < bufferedQueries.size(); i=i+bSize) {
    std::vector<long> tmpBatch;
    for (int j = 0; j < bSize; j++) {
      tmpBatch.push_back(bufferedQueries[i+j]);
    }
    if (shouldDelay) {
      // cout << "with delaying\n";
      // Delayed batching
      vector<int> dist_to_high;
      long total_delays = 0;
      for (int j = 0; j < tmpBatch.size(); j++) {
        // cout << "q" << j << " to highest deg vtx: " << distances[tmpBatch[j]] << endl;
        if (distances[tmpBatch[j]] != MAXLEVEL) {
          dist_to_high.push_back(distances[tmpBatch[j]]);
        } else {
          dist_to_high.push_back(-1);
        }
      }
      int max_dist_to_high = *max_element(dist_to_high.begin(), dist_to_high.end());

      for (int j = 0; j < dist_to_high.size(); j++) {
        if (dist_to_high[j] == -1) {
          dist_to_high[j] = max_dist_to_high;
        }
      }

      for (int j = 0; j < dist_to_high.size(); j++) {
        dist_to_high[j] = max_dist_to_high - dist_to_high[j];
        total_delays += dist_to_high[j];
      }
      // cout << "Total delays (delta): " << total_delays << endl;

      timer t_delay;
      t_delay.start();
      pair<size_t, size_t> share_cnt = Compute_Delay_Skipping(G,tmpBatch,P,dist_to_high);
      t_delay.stop();
      double time1 = t_delay.totalTime;
      batching_time += time1;
      res.push_back(share_cnt);
    } else {
      timer t_t1;
      t_t1.start();
      pair<size_t, size_t> share_cnt = Compute_Base_Skipping(G,tmpBatch,P,true);
      t_t1.stop();
      double time1 = t_t1.totalTime;
      batching_time += time1;
      res.push_back(share_cnt);
    }
  }
  string outstr = "Glign";
  if (selection == 1) outstr = "Sequential";
  if (selection == 2 && !shouldDelay) outstr = "Glign-Intra";  
  if (selection == 2 && shouldDelay) outstr = "Glign-Inter";
  if (selection == 3 && !shouldDelay) outstr = "Glign-Batch";
  if (selection == 3 && shouldDelay) outstr = "Glign";
  cout << outstr + " evaluation time: " << batching_time << endl;
  return res;
}

void ligra_c(int argc, char* argv[]) {
  commandLine P(argc,argv," [-s] <inFile>");
  char* iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");
  bool compressed = P.getOptionValue("-c");
  bool binary = P.getOptionValue("-b");
  bool mmap = P.getOptionValue("-m");
  //cout << "mmap = " << mmap << endl;
  long rounds = P.getOptionLongValue("-rounds",1);
  bool shouldShuffle = P.getOptionValue("-shuffle");

  string queryFileName = string(P.getOptionValue("-qf", ""));
  int combination_max = P.getOptionLongValue("-max_combination", 256);
  size_t bSize = P.getOptionLongValue("-batch", 4);
  size_t n_high_deg = P.getOptionLongValue("-nhighdeg", 4);

  cout << "graph file name: " << iFile << endl;
  cout << "query file name: " << queryFileName << endl;

  // Initialization and preprocessing
  std::vector<long> userQueries; 
  long start = -1;
  char inFileName[300];
  ifstream inFile;
  sprintf(inFileName, queryFileName.c_str());
  inFile.open(inFileName, ios::in);
  while (inFile >> start) {
    userQueries.push_back(start);
  }
  inFile.close();
  // randomly shuffled each run
  std::random_device rd;
  auto rng = std::default_random_engine { rd() };
  // auto rng = std::default_random_engine {};
  if (shouldShuffle)
    std::shuffle(std::begin(userQueries), std::end(userQueries), rng);
  cout << "number of random queries: " << userQueries.size() << endl;

  int setSize = userQueries.size();
  std::vector<long> testQuery[setSize];

  if (symmetric) {
    cout << "symmetric graph\n";
    graph<symmetricVertex> G =
      readGraph<symmetricVertex>(iFile,compressed,symmetric,binary,mmap); //symmetric graph
    // for(int r=0;r<rounds;r++) {
    cout << "n=" << G.n << " m=" << G.m << endl;
    // Streaming...
    vector<long> sortedQueries;
    vector<long> truncatedQueries;
    vector<long> propertySortedQueries;
    tie(truncatedQueries, sortedQueries) = streamingPreprocessing(G, userQueries, n_high_deg, combination_max, P);
    
    // start streaming.
    // input: G, P, bufferedQueries, batch size
    bufferStreaming(G, truncatedQueries, bSize, P);
  } else {
    // For directed graph...
    cout << "asymmetric graph\n";
    graph<asymmetricVertex> G =
      readGraph<asymmetricVertex>(iFile,compressed,symmetric,binary,mmap); //asymmetric graph
    cout << "n=" << G.n << " m=" << G.m << endl;
    
    // Streaming...
    vector<long> sortedQueries;
    vector<long> truncatedQueries;
    vector<long> propertySortedQueries;
    tie(truncatedQueries, sortedQueries) = streamingPreprocessing(G, userQueries, n_high_deg, combination_max, P);
    
    // start streaming.
    // input: G, P, bufferedQueries, batch size
    bufferStreaming(G, truncatedQueries, bSize, P);
  }
}

// for testing CGQ
void test_cgq(int argc, char* argv[]) {
  commandLine P(argc,argv," [-s] <inFile>");
  char* iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");
  bool compressed = P.getOptionValue("-c");
  bool binary = P.getOptionValue("-b");
  bool mmap = P.getOptionValue("-m");
  //cout << "mmap = " << mmap << endl;
  long rounds = P.getOptionLongValue("-rounds",1);
  bool shouldShuffle = P.getOptionValue("-shuffle");


  string queryFileName = string(P.getOptionValue("-qf", ""));
  int combination_max = P.getOptionLongValue("-max_combination", 256);
  size_t bSize = P.getOptionLongValue("-batch", 4);
  size_t n_high_deg = P.getOptionLongValue("-nhighdeg", 4);

  cout << "graph file name: " << iFile << endl;
  cout << "query file name: " << queryFileName << endl;

  // Initialization and preprocessing
  std::vector<long> userQueries; 
  long start = -1;
  char inFileName[300];
  ifstream inFile;
  sprintf(inFileName, queryFileName.c_str());
  inFile.open(inFileName, ios::in);
  while (inFile >> start) {
    userQueries.push_back(start);
  }
  inFile.close();
  // randomly shuffled each run
  std::random_device rd;
  auto rng = std::default_random_engine { rd() };
  // auto rng = std::default_random_engine {};
  if (shouldShuffle)
    std::shuffle(std::begin(userQueries), std::end(userQueries), rng);
  cout << "number of random queries: " << userQueries.size() << endl;

  int setSize = userQueries.size();
  std::vector<long> testQuery[setSize];

  if (symmetric) {
    cout << "symmetric graph\n";
  } else {
    // For directed graph...
    cout << "asymmetric graph\n";
    graph<asymmetricVertex> G =
      readGraph<asymmetricVertex>(iFile,compressed,symmetric,binary,mmap); //asymmetric graph
    cout << "n=" << G.n << " m=" << G.m << endl;
    
    // Streaming...
    vector<long> sortedQueries;
    vector<long> truncatedQueries;
    vector<long> propertySortedQueries;
    tie(truncatedQueries, sortedQueries) = streamingPreprocessing(G, userQueries, n_high_deg, combination_max, P);
    
    // start streaming.
    // input: G, P, bufferedQueries, batch size
    long selection = P.getOptionLongValue("-order",1);
    vector<pair<size_t,size_t>> share1;
    vector<pair<size_t,size_t>> share_unsorted;
    vector<pair<size_t,size_t>> share_sorted;

    set<set<long>> C_Set;
    vector<long> chunk_lookup;

    if (selection == 2) {
      cout << "\non the concurrent setting..\n";
      share_unsorted = streamingWithCGQ(G, truncatedQueries, C_Set, chunk_lookup, bSize, P, true);
    }
  }
}

// skipping status checking
void glign(int argc, char* argv[]) {
  commandLine P(argc,argv," [-s] <inFile>");
  char* iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");  // doesn't affect correctness.
  bool compressed = P.getOptionValue("-c");
  bool binary = P.getOptionValue("-b");
  bool mmap = P.getOptionValue("-m");
  long rounds = P.getOptionLongValue("-rounds",1);  // number of rounds to repeat
  bool shouldShuffle = P.getOptionValue("-shuffle");
  string queryFileName = string(P.getOptionValue("-qf", ""));
  int combination_max = P.getOptionLongValue("-max_combination", 256);
  size_t bSize = P.getOptionLongValue("-batch", 4);
  size_t n_high_deg = P.getOptionLongValue("-nhighdeg", 4);
  long selection = P.getOptionLongValue("-mode",1);
  bool shouldDelay = P.getOptionValue("-delay");

  cout << "graph file name: " << iFile << endl;
  cout << "query file name: " << queryFileName << endl;

  // Initialization and preprocessing
  std::vector<long> userQueries; 
  long start = -1;
  char inFileName[300];
  ifstream inFile;
  sprintf(inFileName, queryFileName.c_str());
  inFile.open(inFileName, ios::in);
  while (inFile >> start) {
    userQueries.push_back(start);
  }
  inFile.close();
  // randomly shuffled each run
  std::random_device rd;
  auto rng = std::default_random_engine { rd() };
  // auto rng = std::default_random_engine {};
  if (shouldShuffle)
    std::shuffle(std::begin(userQueries), std::end(userQueries), rng);
  cout << "number of random queries: " << userQueries.size() << endl;

  int setSize = userQueries.size();
  std::vector<long> testQuery[setSize];

  vector<pair<size_t,size_t>> share1;
  vector<pair<size_t,size_t>> share_unsorted;
  vector<pair<size_t,size_t>> share_sorted;
  vector<size_t> total_N;

  string outstr = "Glign";
  if (selection == 1) outstr = "Sequential";
  if (selection == 2 && !shouldDelay) outstr = "Glign-Intra";  
  if (selection == 2 && shouldDelay) outstr = "Glign-Inter";
  if (selection == 3 && !shouldDelay) outstr = "Glign-Batch";
  if (selection == 3 && shouldDelay) outstr = "Glign";

  if (symmetric) {
    cout << "symmetric graph\n";
    graph<symmetricVertex> G =
      readGraph<symmetricVertex>(iFile,compressed,symmetric,binary,mmap); //symmetric graph
    cout << "n=" << G.n << " m=" << G.m << endl;

    vector<long> sortedQueries;
    vector<long> truncatedQueries;
    vector<long> propertySortedQueries;

    uintE* distances = pbbs::new_array<uintE>(G.n);
    tie(truncatedQueries, sortedQueries) = streamingPreprocessingReturnHops(G, userQueries, n_high_deg, combination_max, P, distances);
    
    // start streaming.
    // input: G, P, bufferedQueries, batch size
    if (selection == 1) {
      cout << "\n" << outstr << " evaluation..\n";
      share1 = bufferStreamingSkipping(G, truncatedQueries, 1, P, distances, true);
      for (int i = 0; i < combination_max; i+=bSize) {
        size_t temp = 0;
        for (int j = i; j < i+bSize; j++) {
          temp += share1[j].first;
        }
        total_N.push_back(temp);
      }
      for (int i = 0; i < total_N.size(); i++) {
        cout << outstr << " F: " << total_N[i] << endl;
      }
    }
    if (selection == 2) {
      cout << "\n" << outstr << " evaluation..\n";
      share_unsorted = bufferStreamingSkipping(G, truncatedQueries, bSize, P, distances, true);
      for (int i = 0; i < share_unsorted.size(); i++) {
        cout << outstr << " F: " << share_unsorted[i].first << endl;
      }
    }
    if (selection == 3) {
      cout << "\n" << outstr << " evaluation..\n";
      share_sorted = bufferStreamingSkipping(G, sortedQueries, bSize, P, distances, true);
      for (int i = 0; i < share_sorted.size(); i++) {
        cout << outstr << " F: " << share_sorted[i].first << endl;
      }
    }

  } else {
    // For directed graph...
    cout << "asymmetric graph\n";
    graph<asymmetricVertex> G =
      readGraph<asymmetricVertex>(iFile,compressed,symmetric,binary,mmap); //asymmetric graph
    cout << "n=" << G.n << " m=" << G.m << endl;
    // Streaming...
    vector<long> sortedQueries;
    vector<long> truncatedQueries;
    vector<long> propertySortedQueries;
    uintE* distances = pbbs::new_array<uintE>(G.n);
    tie(truncatedQueries, sortedQueries) = streamingPreprocessingReturnHops(G, userQueries, n_high_deg, combination_max, P, distances);
    
    // start streaming.
    // input: G, P, bufferedQueries, batch size
    if (selection == 1) {
      cout << "\n" << outstr << " evaluation..\n";
      share1 = bufferStreamingSkipping(G, truncatedQueries, 1, P, distances, true);
      for (int i = 0; i < combination_max; i+=bSize) {
      size_t temp = 0;
      for (int j = i; j < i+bSize; j++) {
        temp += share1[j].first;
      }
      total_N.push_back(temp);
      }
      for (int i = 0; i < total_N.size(); i++) {
        cout << outstr << " F: " <<  total_N[i] << endl;
      }
    }
    if (selection == 2) {
      cout << "\n" << outstr << " evaluation..\n";
      share_unsorted = bufferStreamingSkipping(G, truncatedQueries, bSize, P, distances, true);
      for (int i = 0; i < share_unsorted.size(); i++) {
        cout << outstr << " F: " <<  share_unsorted[i].first << endl;
      }
    }
    if (selection == 3) {
      cout << "\n" << outstr << " evaluation..\n";
      share_sorted = bufferStreamingSkipping(G, sortedQueries, bSize, P, distances, true);
      for (int i = 0; i < share_sorted.size(); i++) {
        cout << outstr << " F: " <<  share_sorted[i].first << endl;
      }
    }
  }
}

// for testing iBFS heuristic
void iBFS(int argc, char* argv[]) {
  commandLine P(argc,argv," [-s] <inFile>");
  char* iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");
  bool compressed = P.getOptionValue("-c");
  bool binary = P.getOptionValue("-b");
  bool mmap = P.getOptionValue("-m");
  //cout << "mmap = " << mmap << endl;
  long rounds = P.getOptionLongValue("-rounds",1);

  string queryFileName = string(P.getOptionValue("-qf", ""));
  int combination_max = P.getOptionLongValue("-max_combination", 512);
  size_t bSize = P.getOptionLongValue("-batch", 4);
  size_t n_high_deg = P.getOptionLongValue("-nhighdeg", 4);

  cout << "graph file name: " << iFile << endl;
  cout << "query file name: " << queryFileName << endl;

  // Initialization and preprocessing
  std::vector<long> userQueries; 
  long start = -1;
  char inFileName[300];
  ifstream inFile;
  sprintf(inFileName, queryFileName.c_str());
  inFile.open(inFileName, ios::in);
  while (inFile >> start) {
    userQueries.push_back(start);
  }
  inFile.close();
  // randomly shuffled each run
  std::random_device rd;
  auto rng = std::default_random_engine { rd() };
  // std::shuffle(std::begin(userQueries), std::end(userQueries), rng);
  cout << "number of random queries: " << userQueries.size() << endl;

  int setSize = userQueries.size();
  std::vector<long> testQuery[setSize];

  if (symmetric) {
    cout << "No need to implement.\n";
  } else {
    // For directed graph...
    cout << "asymmetric graph\n";
    graph<asymmetricVertex> G =
      readGraph<asymmetricVertex>(iFile,compressed,symmetric,binary,mmap); //asymmetric graph
    cout << "n=" << G.n << " m=" << G.m << endl;
    
    // Streaming...
    vector<long> sortedQueries;
    vector<long> truncatedQueries;
    vector<long> propertySortedQueries;

    uintE* distances = pbbs::new_array<uintE>(G.n);

    tie(truncatedQueries, sortedQueries) = streamingPreprocessingReturnHops(G, userQueries, n_high_deg, combination_max, P, distances);
    
    // start streaming.
    // input: G, P, bufferedQueries, batch size
    long selection = P.getOptionLongValue("-mode",1);

    // testing to see if queries share neighbors
    cout << "\nTesting " << truncatedQueries.size() << " queries" << endl;
    int p_seq[4] = {4, 16, 64, 128};
    vector<long> ibfs_sorted;
    set<long> selected_q;
    int vtx_share_cnt[truncatedQueries.size()] = {}; // not needed
    set<long> q_set;
    for (long i = 0; i < truncatedQueries.size(); i++) {
      q_set.insert(truncatedQueries[i]);
    }

    for (int s = 0; s < 4; s++) {
      int p = p_seq[s];
      for (long i = 0; i < truncatedQueries.size(); i++) {
        vtx_share_cnt[i] = 0;
        long _vtx = truncatedQueries[i];
        intE* nghs = G.V[_vtx].getOutNeighbors();
        long out_deg_i = G.V[_vtx].getOutDegree();
        set<long> i_set;

        for (long j = 0; j < out_deg_i; j++) {
          long _ngh = nghs[j];
          i_set.insert(_ngh);
        }
        for (long k = i; k < truncatedQueries.size(); k++) {
          if (k == i) continue;
          long _vtx_k = truncatedQueries[k];
          intE* nghs_k = G.V[_vtx_k].getOutNeighbors();
          long out_deg_k = G.V[_vtx_k].getOutDegree();
          for (long j = 0; j < out_deg_k; j++) {
            long _ngh = nghs_k[j];
            if (i_set.find(_ngh) != i_set.end()) {
              if (G.V[_ngh].getOutDegree() > 128) {
                if (out_deg_i < p && out_deg_k < p) {
                  // cout << i << " and " << k << " share at least one neighbor that q > 128, p < " << p << ", degs: " << out_deg_i << " " << out_deg_k << endl;
                  if (selected_q.find(_vtx) == selected_q.end()) {
                    // cout << i << " share at least one neighbor that q > 128, p < " << p << ", degs: " << out_deg_i << endl;
                    ibfs_sorted.push_back(_vtx);
                    selected_q.insert(_vtx);
                  }
                  if (selected_q.find(_vtx_k) == selected_q.end()) {
                    // cout << k << " share at least one neighbor that q > 128, p < " << p << ", degs: " << out_deg_k << endl;
                    ibfs_sorted.push_back(_vtx_k);
                    selected_q.insert(_vtx_k);
                  }
                }
                vtx_share_cnt[i]++;
                break;
              }
            }
          }
        }
      }
    }
    vector<long> rest_queries;
    vector<long> rest_queries_original;
    vector<long> dummy_queries;
    vector<long> combined_sorted;
    for (long i = 0; i < ibfs_sorted.size(); i++) {
      combined_sorted.push_back(ibfs_sorted[i]);
    }
    for (long i = 0; i < truncatedQueries.size(); i++) {
      if (selected_q.find(truncatedQueries[i]) == selected_q.end()) {
        rest_queries_original.push_back(truncatedQueries[i]);
      }
    }
    tie(dummy_queries, rest_queries) = streamingPreprocessing(G, rest_queries_original, n_high_deg, combination_max-ibfs_sorted.size(), P);
    for (long i = 0; i < rest_queries.size(); i++) {
      combined_sorted.push_back(rest_queries[i]);
    }

    for (long i = 0; i < truncatedQueries.size(); i++) {
      if (selected_q.find(truncatedQueries[i]) == selected_q.end()) {
        ibfs_sorted.push_back(truncatedQueries[i]);
      }
    }

    cout << "==== ibfs ====\n";
    bufferStreamingSkipping(G, ibfs_sorted, bSize, P, distances, true);
  }
}

vector<vector<long>> my_comb(int N, int K, int max_size=LONG_MAX)
{
  std::string bitmask(K, 1); // K leading 1's
  bitmask.resize(N, 0); // N-K trailing 0's
  vector<vector<long>> ret;
  // print integers and permute bitmask
  do {
    vector<long> tmp_sample;
    for (int i = 0; i < N; ++i) // [0..N-1] integers
    {
        if (bitmask[i]) {
            tmp_sample.push_back(i);
        }
    }
    ret.push_back(tmp_sample);
    if (ret.size() >= max_size) {   // break early
        return ret;
    }
  } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
  return ret;
}

int query_generation_skew(int argc, char* argv[]) {
  commandLine P(argc,argv," [-s] <inFile>");
  char* iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");
  bool compressed = P.getOptionValue("-c");
  bool binary = P.getOptionValue("-b");
  bool mmap = P.getOptionValue("-m");
  //cout << "mmap = " << mmap << endl;
  long rounds = P.getOptionLongValue("-rounds",1);
  // long total_len = P.getOptionLongValue("-total",512);
  bool shouldShuffle = P.getOptionValue("-shuffle");

  string queryFileName = string(P.getOptionValue("-qf", ""));
  int combination_max = P.getOptionLongValue("-max_combination", 256);
  size_t bSize = P.getOptionLongValue("-batch", 4);
  size_t n_high_deg = P.getOptionLongValue("-nhighdeg", 4);

  cout << "graph file name: " << iFile << endl;
  cout << "query file name: " << queryFileName << endl;

  // Initialization and preprocessing
  std::vector<long> userQueries; 
  long start = -1;
  char inFileName[300];
  ifstream inFile;
  sprintf(inFileName, queryFileName.c_str());
  inFile.open(inFileName, ios::in);
  while (inFile >> start) {
    userQueries.push_back(start);
  }
  inFile.close();

  // randomly shuffled each run
  std::random_device rd;
  auto rng = std::default_random_engine { rd() };
  // auto rng = std::default_random_engine {};
  // if (shouldShuffle)
    std::shuffle(std::begin(userQueries), std::end(userQueries), rng);

  cout << "number of queries: " << userQueries.size() << endl;
  int batch_size = userQueries.size();

  cout << "asymmetric graph\n";
  graph<asymmetricVertex> G =
    readGraph<asymmetricVertex>(iFile,compressed,symmetric,binary,mmap); //asymmetric graph
  cout << "n=" << G.n << " m=" << G.m << endl;

  vector<long> sortedQueries;
  vector<long> truncatedQueries;
  size_t n = G.n;
  std::vector<std::pair<long, long>> vIDDegreePairs;
  for (long i = 0; i < n; i++) {
    long temp_degree =  G.V[i].getOutDegree();
    if (temp_degree >= 50) {
      vIDDegreePairs.push_back(std::make_pair(i, temp_degree));
    }
  }
  std::sort(vIDDegreePairs.begin(), vIDDegreePairs.end(), sortByLargerSecondElement);
  vector<long> highdegQ;
  int high_deg_batch = n_high_deg;
  for (int i = 0; i < high_deg_batch; i++) {
    highdegQ.push_back(vIDDegreePairs[i].first);
  }
  uintE* distances = pbbs::new_array<uintE>(n);
  uintE* distances_multiple;
  // On edge reversed graph...
  G.transpose();
  distances_multiple = Compute_Eval(G,highdegQ,P);  // to get hops (bfs distance)
  // uintE* distances = pbbs::new_array<uintE>(n);
  G.transpose();
  parallel_for(size_t i = 0; i < n; i++) {
    distances[i] = (uintE)MAXLEVEL;
  }
  parallel_for(size_t i = 0; i < n; i++) {
    for (int j = 0; j < high_deg_batch; j++) {
      if (distances_multiple[j+i*high_deg_batch] < distances[i]) {
        distances[i] = distances_multiple[j+i*high_deg_batch];
      }
    }
  }

  // hop distributions of input queries.
  std::map<long, long> user_hist;
  for (long i = 0; i < userQueries.size(); i++) {
    int dist = distances[userQueries[i]];
    user_hist[dist]++;
  }
  for (const auto& x : user_hist) std::cout << x.first << " " << x.second <<"\n";

  std::vector<std::pair<size_t, long>> vtxHopPairs;
  cout << "\ninput queries: \n";
  for (long i = 0; i < userQueries.size(); i++) {
    // cout << truncatedQueries[i] << ": " << distances[truncatedQueries[i]] << endl;
    vtxHopPairs.push_back(std::make_pair(userQueries[i], distances[userQueries[i]]));
  }
  std::sort(vtxHopPairs.begin(), vtxHopPairs.end(), sortByLargerSecondElement);

  int _count = 0;
  for (int i = vtxHopPairs.size()-1; i >=0; i--) {
    if (vtxHopPairs[i].second == 1) {
      if (_count < 96) {
        sortedQueries.push_back(vtxHopPairs[i].first);
        _count++;
      } else {
        break;
      }
    }
  }
  _count = 0;
  int last = 0;
  for (int i = vtxHopPairs.size()-1; i >=0; i--) {
    if (vtxHopPairs[i].second == 2) {
      if (_count < 96) {
        sortedQueries.push_back(vtxHopPairs[i].first);
        _count++;
        last = i;
      } else {
        break;
      }
    }
  }
  cout << "last: " << last << endl;

  string outFileName = string(P.getOptionValue("-of", ""));
  ofstream outputFile (outFileName, ios::out | ios::binary);
  if (!outputFile.is_open()) {
    std::cout << "Unable to open file: " << outFileName << std::endl;
    return -1;
  }
  srand((unsigned)time(NULL));
  for (int k = 0; k < combination_max - 192; k++) {
    size_t vidx = rand() % last;
    // cout << vtxHopPairs[vidx].first << ": " << vtxHopPairs[vidx].second << endl;
    sortedQueries.push_back(vtxHopPairs[vidx].first);
  }

  cout << "\nskewed queries: \n";
  std::shuffle(std::begin(sortedQueries), std::end(sortedQueries), rng);
  for (long i = 0; i < sortedQueries.size(); i++) {
    cout << sortedQueries[i] << ": " << distances[sortedQueries[i]] << endl;
    // cout << sortedQueries[i] << endl;
    outputFile << sortedQueries[i] << endl;
  }

  outputFile.close();
  return 0;
}

int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv," [-s] <inFile>");
  string options = string(P.getOptionValue("-option", "scenario3"));
    
  if (options == "ligra-c") {
    cout << "running Ligra-C\n";
    ligra_c(argc, argv);  // reordering
  }

  if (options == "ibfs") {
    cout << "running iBFS heuristic\n";
    iBFS(argc, argv);
  }

  if (options == "ground-truth") {
    cout << "Table 13: Ground Truth Study\n";
    groundtruth(argc, argv);
  }

  if (options == "glign") {
    cout << "running batched queries on glign\n";
    glign(argc, argv);
  }

  // if (options == "cgq") {
  //   cout << "testing async concurrent query processing on Ligra.\n";
  //   test_cgq(argc, argv);
  // }
}
#endif
