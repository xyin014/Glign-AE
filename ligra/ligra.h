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
#define MAXPATH 10000
#define MAXLEVEL 10000
#define MAXWIDTH 10000

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
    cout << "dense mod\n";
    if(degrees) free(degrees);
    if(frontierVertices) free(frontierVertices);
    vs.toDense();
    return (fl & dense_forward) ?
      edgeMapDenseForward<data, vertex, VS, F>(GA, vs, f, fl) :
      edgeMapDense<data, vertex, VS, F>(GA, vs, f, fl);
  } else {
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
template<class vertex>
uintE* Compute_Eval_Prop(graph<vertex>&, vector<long>, commandLine);
template<class vertex>
void Compute_Base(graph<vertex>&, vector<long>, commandLine, bool should_profile=false);
template<class vertex>
void Compute_Delay(graph<vertex>&, vector<long>, commandLine, vector<int>, bool should_profile=false);

template<class vertex>
void Compute(hypergraph<vertex>&, commandLine);

bool sortByLargerSecondElement(const pair<long, long> &a, const pair<long, long> &b) {
  return (a.second > b.second);
}

// for query generation
int QueryGeneration_Random(int argc, char* argv[]) {
  commandLine P(argc,argv," [-s] <inFile>");
  char* iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");
  bool compressed = P.getOptionValue("-c");
  bool binary = P.getOptionValue("-b");
  bool mmap = P.getOptionValue("-m");

  string outFileName = string(P.getOptionValue("-of", ""));
  ofstream outputFile (outFileName, ios::out | ios::binary);
  if (!outputFile.is_open()) {
    std::cout << "Unable to open file: " << outFileName << std::endl;
    return -1;
  }

  size_t n_high_deg = P.getOptionLongValue("-nhighdeg", 4);
  int q_size = P.getOptionIntValue("-querysize", 8192*4);
  std::vector<size_t> queries;
  unordered_set <size_t> querySet;
  srand((unsigned)time(NULL));

  if (symmetric) {
    cout << "symmetric graph\n";
    graph<symmetricVertex> G =
      readGraph<symmetricVertex>(iFile,compressed,symmetric,binary,mmap); //symmetric graph
    // for(int r=0;r<rounds;r++) {
    cout << "n=" << G.n << " m=" << G.m << endl;
    size_t n = G.n;

    while (querySet.size() < q_size) {
      size_t vtx = rand() % n;
      if (querySet.find(vtx) == querySet.end()) {
        querySet.insert(vtx);
        cout << vtx <<  endl;
        outputFile << vtx << endl;
      }
    }
    outputFile.close();
  } else {
    cout << "asymmetric graph\n";
    graph<asymmetricVertex> G =
      readGraph<asymmetricVertex>(iFile,compressed,symmetric,binary,mmap); //asymmetric graph
    cout << "n=" << G.n << " m=" << G.m << endl;
    size_t n = G.n;

    while (querySet.size() < q_size) {
      size_t vtx = rand() % n;
      if (querySet.find(vtx) == querySet.end()) {
        querySet.insert(vtx);
        cout << vtx <<  endl;
        outputFile << vtx << endl;
      }
    }
    outputFile.close();
  }

  return 0;
}

// for query generation
int QueryGeneration_Hops(int argc, char* argv[]) {

  return 0;
}

// for query generation
int QueryGeneration_Property(int argc, char* argv[]) {
  commandLine P(argc,argv," [-s] <inFile>");
  char* iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");
  bool compressed = P.getOptionValue("-c");
  bool binary = P.getOptionValue("-b");
  bool mmap = P.getOptionValue("-m");
  
  string outFileName = string(P.getOptionValue("-of", ""));
  ofstream outputFile (outFileName, ios::out | ios::binary);
  if (!outputFile.is_open()) {
    std::cout << "Unable to open file: " << outFileName << std::endl;
    return -1;
  }

  size_t n_high_deg = P.getOptionLongValue("-nhighdeg", 4);
  int q_size = P.getOptionIntValue("-querysize", 8192*4);
  std::vector<size_t> queries;
  unordered_set <size_t> querySet;
  srand((unsigned)time(NULL));

  // if (symmetric) {
  //   cout << "symmetric graph\n";
  //   graph<symmetricVertex> G =
  //     readGraph<symmetricVertex>(iFile,compressed,symmetric,binary,mmap); //symmetric graph
  //   // for(int r=0;r<rounds;r++) {
  //   cout << "n=" << G.n << " m=" << G.m << endl;
  //   size_t n = G.n;

  //   // finding high degree vtxs.
  //   std::vector<std::pair<long, long>> vIDDegreePairs;
  //   for (long i = 0; i < n; i++) {
  //     long temp_degree =  G.V[i].getOutDegree();
  //     if (temp_degree >= 50) {
  //       vIDDegreePairs.push_back(std::make_pair(i, temp_degree));
  //     }
  //   }
  //   std::sort(vIDDegreePairs.begin(), vIDDegreePairs.end(), sortByLargerSecondElement);
  //   vector<long> highdegQ;
  //   int high_deg_batch = n_high_deg;
  //   for (int i = 0; i < high_deg_batch; i++) {
  //     highdegQ.push_back(vIDDegreePairs[i].first);
  //   }

  //   uintE* distances_multiple;
  //   distances_multiple = Compute_Eval_Prop(G,highdegQ,P);
  //   uintE* distances = pbbs::new_array<uintE>(n);
  //   parallel_for(size_t i = 0; i < n; i++) {
  //     distances[i] = (uintE)MAXLEVEL;
  //   }
  //   parallel_for(size_t i = 0; i < n; i++) {
  //     for (int j = 0; j < high_deg_batch; j++) {
  //       if (distances_multiple[j+i*high_deg_batch] < distances[i]) {
  //         distances[i] = distances_multiple[j+i*high_deg_batch];
  //       }
  //     }
  //   }


    
  //   outputFile.close();
  // } else {
  //   cout << "asymmetric graph\n";
  //   graph<asymmetricVertex> G =
  //     readGraph<asymmetricVertex>(iFile,compressed,symmetric,binary,mmap); //asymmetric graph
  //   cout << "n=" << G.n << " m=" << G.m << endl;
  //   size_t n = G.n;

    
  //   outputFile.close();
  // }

  return 0;
}

template <class vertex>
void bufferStreaming_with_arrival_time(graph<vertex>& G, std::vector<long> bufferedQueries, std::vector<double> arrivalTimes, uintE* distances, int bSize, commandLine P, bool should_profile=false) {
  vector<double> latency_map;

  double static_latency = 0.0;
  // double start_time1 = update_timer.get_time();
  double earlier_start1 = 0;
  double new_est = 0;
  // double new_est = arrivalTimes[0+bSize-1];
  timer start_time1; start_time1.start();
  for (int i = 0; i < bufferedQueries.size(); i=i+bSize) {
    cout << "i: " << i << endl;
    std::vector<long> tmpBatch;
    double arrival_last_in_the_batch = arrivalTimes[i+bSize-1]; // last arrived in the batch
    for (int j = 0; j < bSize; j++) {
      tmpBatch.push_back(bufferedQueries[i+j]);
      if (arrival_last_in_the_batch < new_est) {
        // cout << "new_est - arrivalTimes[i+j]: " << new_est - arrivalTimes[i+j] << endl;
        static_latency += new_est - arrivalTimes[i+j];
      } else {
        static_latency += arrival_last_in_the_batch - arrivalTimes[i+j];
      }
    }

    // for delaying
    vector<int> dist_to_high;
    long total_delays = 0;
    for (int j = 0; j < tmpBatch.size(); j++) {
      cout << "q" << j << " to highest deg vtx: " << distances[tmpBatch[j]] << endl;
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
      cout << "No. " << j << " defer " << dist_to_high[j] << " iterations\n";
    }
    cout << "Total delays (delta): " << total_delays << endl;

    timer t_t1;
    t_t1.start();
    Compute_Delay(G,tmpBatch,P,dist_to_high);
    t_t1.stop();
    double time1 = t_t1.totalTime;

    // record latency for each query
    for (int ii = 0; ii < bSize; ii++) {
      latency_map.push_back(new_est+time1);
    }

    if (arrival_last_in_the_batch < new_est) {
      new_est = time1 + new_est;
    } else {
      new_est = time1 + new_est + arrival_last_in_the_batch - new_est;
    }
    static_latency += (time1)*bSize;
    
    cout << "current latency: " << static_latency << endl;
  }
  start_time1.stop();
  double query_time1 = start_time1.totalTime;
  cout << "static batching version query time: " << query_time1 << endl;
  cout << "Static total latency: " << static_latency << endl;

  double check_sum = 0.0;
  sort(latency_map.begin(), latency_map.end());
  for (int ii = 0; ii < bufferedQueries.size(); ii++) {
    cout << latency_map[ii] << endl;
    check_sum += latency_map[ii];
  }
  cout << "check_sum: " << check_sum << endl;

  // if (should_profile) {
  //   for (int i = 0; i < bufferedQueries.size(); i=i+bSize) {
  //     std::vector<long> tmpBatch;
  //     for (int j = 0; j < bSize; j++) {
  //       tmpBatch.push_back(bufferedQueries[i+j]);
  //     }
  //     Compute_Base(G,tmpBatch,P,true);
  //   }
  // }
}

void scenario_adaptive(int argc, char* argv[]) {
  commandLine P(argc,argv," [-s] <inFile>");
  char* iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");
  bool compressed = P.getOptionValue("-c");
  bool binary = P.getOptionValue("-b");
  bool mmap = P.getOptionValue("-m");
  long rounds = P.getOptionLongValue("-rounds",1);

  string queryFileName = string(P.getOptionValue("-qf", ""));
  string arrivalFileName = string(P.getOptionValue("-af", ""));
  int combination_max = P.getOptionLongValue("-max_combination", 256);
  size_t bSize = P.getOptionLongValue("-batch", 4);
  size_t n_high_deg = P.getOptionLongValue("-nhighdeg", 4);

  cout << "graph file name: " << iFile << endl;
  cout << "query file name: " << queryFileName << endl;

  // Initialization and preprocessing
  std::vector<long> userQueries;
  std::vector<double> arrivalTimes;
  long start = -1;
  char inFileName[300];
  char inFileName2[300];
  ifstream inFile;
  sprintf(inFileName, queryFileName.c_str());
  inFile.open(inFileName, ios::in);
  while (inFile >> start) {
    userQueries.push_back(start);
  }
  inFile.close();

  ifstream inFile2;
  sprintf(inFileName2, arrivalFileName.c_str());
  inFile2.open(inFileName2, ios::in);
  double aa, ff;
  while (inFile2 >> aa >> ff) {
    arrivalTimes.push_back(aa);
  }
  inFile2.close();

  // randomly shuffled each run
  std::random_device rd;
  auto rng = std::default_random_engine { rd() };
  // auto rng = std::default_random_engine {};
  // std::shuffle(std::begin(userQueries), std::end(userQueries), rng);
  cout << "number of random queries: " << userQueries.size() << endl;

  int setSize = userQueries.size();
  std::vector<long> testQuery[setSize];

  if (symmetric) {
    cout << "symmetric graph\n";
    graph<symmetricVertex> G =
      readGraph<symmetricVertex>(iFile,compressed,symmetric,binary,mmap); //symmetric graph
    // for(int r=0;r<rounds;r++) {
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
    for (int i = 0; i < high_deg_batch; i++) {
      highdegQ.push_back(vIDDegreePairs[i].first);
    }

    uintE* distances_multiple;
    distances_multiple = Compute_Eval(G,highdegQ,P);
    uintE* distances = pbbs::new_array<uintE>(n);
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

    // start streaming.
    // input: G, P, bufferedQueries, batch size
    cout << "\non the unsorted buffer..\n";
    bufferStreaming_with_arrival_time(G, userQueries, arrivalTimes, distances, bSize, P);
    cout << endl;

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
    for (int i = 0; i < high_deg_batch; i++) {
      highdegQ.push_back(vIDDegreePairs[i].first);
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

    // Streaming...
    cout << "\non the unsorted buffer..\n";
    bufferStreaming_with_arrival_time(G, userQueries, arrivalTimes, distances, bSize, P);
    cout << endl;

  }

}

// for dealying
void scenario3(int argc, char* argv[]) {
  commandLine P(argc,argv," [-s] <inFile>");
  char* iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");
  bool compressed = P.getOptionValue("-c");
  bool binary = P.getOptionValue("-b");
  bool mmap = P.getOptionValue("-m");
  //cout << "mmap = " << mmap << endl;
  long rounds = P.getOptionLongValue("-rounds",1);

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
  std::shuffle(std::begin(userQueries), std::end(userQueries), rng);

  cout << "number of random queries: " << userQueries.size() << endl;
  int batch_size = userQueries.size();
  
  if (symmetric) {
    cout << "symmetric graph\n";
    graph<symmetricVertex> G =
      readGraph<symmetricVertex>(iFile,compressed,symmetric,binary,mmap); //symmetric graph
    // for(int r=0;r<rounds;r++) {
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
    distances_multiple = Compute_Eval(G,highdegQ,P);
    uintE* distances = pbbs::new_array<uintE>(n);
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
    for (int i = 0; i < combination_max; i++) {
      timer t_seq, t_batch, t_delay;
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
      t_seq.start();
      for (int j = 0; j < tmp_batch.size(); j++) {
        vector<long> tmp_single_query;
        tmp_single_query.push_back(tmp_batch[j]);
        Compute_Base(G,tmp_single_query,P);
      }
      t_seq.stop();

      // Batching
      t_batch.start();
      Compute_Base(G,tmp_batch,P);
      t_batch.stop();


      // Delayed batching
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
        cout << "No. " << j << " defer " << dist_to_high[j] << " iterations\n";
      }
      cout << "Total delays (delta): " << total_delays << endl;

      t_delay.start();
      Compute_Delay(G,tmp_batch,P,dist_to_high);
      t_delay.stop();

      double seq_time = t_seq.totalTime;
      double batch_time = t_batch.totalTime;
      double delay_time = t_delay.totalTime;
      t_seq.reportTotal("sequential time");
      t_batch.reportTotal("batching evaluation time");
      t_delay.reportTotal("delayed batching evaluation time");

      cout << "Batching speedup: " << seq_time / batch_time << endl;
      cout << "Delayed batching speedup: " << seq_time / delay_time << endl;

      // profiling the affinities
      Compute_Base(G,tmp_batch,P,true);
      Compute_Delay(G,tmp_batch,P,dist_to_high,true);

      cout << "=================\n";
    }
      
    // }
    G.del();
    pbbs::delete_array(distances, n);
    pbbs::delete_array(distances_multiple, n*highdegQ.size());

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
    for (int i = 0; i < combination_max; i++) {
      timer t_seq, t_batch, t_delay;
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
      t_seq.start();
      for (int j = 0; j < tmp_batch.size(); j++) {
        vector<long> tmp_single_query;
        tmp_single_query.push_back(tmp_batch[j]);
        Compute_Base(G,tmp_single_query,P);
      }
      t_seq.stop();

      // Batching
      vector<int> dummy_dist;
      for (int j = 0; j < tmp_batch.size(); j++) {
        dummy_dist.push_back(0);
      }
      t_batch.start();
      for (int r = 0; r < rounds; r++) {
        Compute_Base(G,tmp_batch,P);
        // Compute_Delay(G,tmp_batch,P,dummy_dist);
      }
      t_batch.stop();

      // Delayed batching
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
        cout << "No. " << j << " defer " << dist_to_high[j] << " iterations\n";
      }
      cout << "Total delays (delta): " << total_delays << endl;

      t_delay.start();
      for (int r = 0; r < rounds; r++) {
        Compute_Delay(G,tmp_batch,P,dist_to_high);
        // Compute_Delay(G,tmp_batch,P,dummy_dist);
        // Compute_Base(G,tmp_batch,P);
      }
      t_delay.stop();

      // // Batching
      // t_batch.start();
      // for (int r = 0; r < rounds; r++) {
      //   Compute_Base(G,tmp_batch,P);
      // }
      // t_batch.stop();

      // // Sequential
      // t_seq.start();
      // for (int j = 0; j < tmp_batch.size(); j++) {
      //   vector<long> tmp_single_query;
      //   tmp_single_query.push_back(tmp_batch[j]);
      //   Compute_Base(G,tmp_single_query,P);
      // }
      // t_seq.stop();

      double seq_time = t_seq.totalTime;
      double batch_time = t_batch.totalTime / rounds;
      double delay_time = t_delay.totalTime / rounds;
      t_seq.reportTotal("sequential time");
      t_batch.reportTotal("batching evaluation time");
      t_delay.reportTotal("delayed batching evaluation time");

      cout << "Batching speedup: " << seq_time / batch_time << endl;
      cout << "Delayed batching speedup: " << seq_time / delay_time << endl;

      // profiling the affinities
      Compute_Base(G,tmp_batch,P,true);
      Compute_Delay(G,tmp_batch,P,dist_to_high,true);

      // // To do: remove
      // for (int j = 0; j < tmp_batch.size(); j++) {
      //   vector<long> tmp_single_query;
      //   tmp_single_query.push_back(tmp_batch[j]);
      //   Compute_Base(G,tmp_single_query,P, true);
      // }
      cout << "=================\n";
    }
      
    // }
    G.del();
    pbbs::delete_array(distances, n);
    pbbs::delete_array(distances_multiple, n*highdegQ.size());
  }
}

void scenario_adaptive_new(int argc, char* argv[]) {
  commandLine P(argc,argv," [-s] <inFile>");
  char* iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");
  bool compressed = P.getOptionValue("-c");
  bool binary = P.getOptionValue("-b");
  bool mmap = P.getOptionValue("-m");
  //cout << "mmap = " << mmap << endl;
  long rounds = P.getOptionLongValue("-rounds",1);
  bool seq = P.getOptionValue("-i");

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

  cout << "number of random queries: " << userQueries.size() << endl;
  int batch_size = userQueries.size();
  
  if (symmetric) {
    cout << "symmetric graph\n";
    graph<symmetricVertex> G =
      readGraph<symmetricVertex>(iFile,compressed,symmetric,binary,mmap); //symmetric graph
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
    distances_multiple = Compute_Eval(G,highdegQ,P);
    uintE* distances = pbbs::new_array<uintE>(n);
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
    for (int i = 0; i < batchedQuery.size(); i=i+bSize) {
      timer t_seq, t_batch, t_delay;
      vector<long> tmp_batch;
      cout << "Evaluating queries: ";
      for (int j = 0; j < bSize; j++) {
        long tmp_query_id = batchedQuery[j];
        tmp_batch.push_back(tmp_query_id);
        cout << tmp_query_id << " ";
      }
      cout << endl;

      // Sequential
      // t_seq.start();
      // for (int j = 0; j < tmp_batch.size(); j++) {
      //   vector<long> tmp_single_query;
      //   tmp_single_query.push_back(tmp_batch[j]);
      //   Compute_Base(G,tmp_single_query,P);
      // }
      // t_seq.stop();

      // Delayed batching
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
        cout << "No. " << j << " defer " << dist_to_high[j] << " iterations\n";
      }
      cout << "Total delays (delta): " << total_delays << endl;

      t_delay.start();
      Compute_Delay(G,tmp_batch,P,dist_to_high);
      t_delay.stop();

      // double seq_time = t_seq.totalTime;
      double delay_time = t_delay.totalTime;
      // t_seq.reportTotal("sequential time");
      t_delay.reportTotal("delayed batching evaluation time");

      cout << "=================\n";
    }
    G.del();
    pbbs::delete_array(distances, n);
    pbbs::delete_array(distances_multiple, n*highdegQ.size());

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
    for (int i = 0; i < batchedQuery.size(); i=i+bSize) {
      timer t_seq, t_batch, t_delay;
      vector<long> tmp_batch;
      cout << "Evaluating queries: ";
      for (int j = 0; j < bSize; j++) {
        long tmp_query_id = batchedQuery[i+j];
        tmp_batch.push_back(tmp_query_id);
        cout << tmp_query_id << " ";
      }
      cout << endl;

      // // Sequential
      // t_seq.start();
      // for (int j = 0; j < tmp_batch.size(); j++) {
      //   vector<long> tmp_single_query;
      //   tmp_single_query.push_back(tmp_batch[j]);
      //   Compute_Base(G,tmp_single_query,P);
      // }
      // t_seq.stop();

      // Delayed batching
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
        cout << "No. " << j << " defer " << dist_to_high[j] << " iterations\n";
      }
      cout << "Total delays (delta): " << total_delays << endl;

      t_delay.start();
      for (int r = 0; r < rounds; r++) {
        Compute_Delay(G,tmp_batch,P,dist_to_high);
      }
      t_delay.stop();


      // double seq_time = t_seq.totalTime;
      double delay_time = t_delay.totalTime / rounds;
      // t_seq.reportTotal("sequential time");
      t_delay.reportTotal("delayed batching evaluation time");

      cout << "=================\n";
    }
      
    // }
    G.del();
    pbbs::delete_array(distances, n);
    pbbs::delete_array(distances_multiple, n*highdegQ.size());
  }
}

void test_4(int argc, char* argv[]) {
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

  // randomly shuffled each run
  std::random_device rd;
  auto rng = std::default_random_engine { rd() };
  // auto rng = std::default_random_engine {};
  std::shuffle(std::begin(userQueries), std::end(userQueries), rng);

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
      t_seq.start();
      for (int j = 0; j < tmp_batch.size(); j++) {
        vector<long> tmp_single_query;
        tmp_single_query.push_back(tmp_batch[j]);
        Compute_Base(G,tmp_single_query,P);
      }
      t_seq.stop();
      double seq_time = t_seq.totalTime;

      // Batching
      t_batch.start();
      for (int r = 0; r < rounds; r++) {
        Compute_Base(G,tmp_batch,P);
        // Compute_Delay(G,tmp_batch,P,dummy_dist);
      }
      t_batch.stop();

      vector<int> tmp_delay = {0, 0};
      int max_ind = distances[tmp_batch[0]] > distances[tmp_batch[1]] ? 0 : 1;
      int min_ind = distances[tmp_batch[0]] < distances[tmp_batch[1]] ? 0 : 1;
      vector<double> combined_spd;
      for (int j = 0; j < max_delay; j++) {
        tmp_delay[min_ind] = j;
        cout << "No. " << min_ind << " delays " << j << " iterations\n";
        timer t_delay; t_delay.start();
        Compute_Delay(G,tmp_batch,P,tmp_delay);
        // Compute_Base(G,tmp_batch,P);
        t_delay.stop();
        double delay_time = t_delay.totalTime;
        combined_spd.push_back(seq_time / delay_time);
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
        cout << "No. " << j << " defer " << dist_to_high[j] << " iterations\n";
      }
      cout << "Total delays (delta): " << total_delays << endl;

      // t_delay.start();
      // for (int r = 0; r < rounds; r++) {
      //   Compute_Delay(G,tmp_batch,P,dist_to_high);
      //   // Compute_Base(G,tmp_batch,P);
      // }
      // t_delay.stop();


      
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

      int best_delay = -1;
      double best_spd = 0.0;
      for (int j = 0; j < combined_spd.size(); j++) {
        if (combined_spd[j] > best_spd) {
          best_spd = combined_spd[j];
          best_delay = j;
        }
      }
      int tmp_delta = abs(dist_to_high[0] - dist_to_high[1]);
      cout << "Best delay: " << best_delay << " " << best_spd  << endl;
      cout << "Delta_Speedup: " << tmp_delta << " " << combined_spd[tmp_delta] << endl;


      // To do: remove
      for (int j = 0; j < tmp_batch.size(); j++) {
        vector<long> tmp_single_query;
        tmp_single_query.push_back(tmp_batch[j]);
        Compute_Base(G,tmp_single_query,P, true);
      }

      cout << "=================\n";
    }
      
    // }
    G.del();
    pbbs::delete_array(distances, n);
    pbbs::delete_array(distances_multiple, n*highdegQ.size());
  }
}

void test_5(int argc, char* argv[]) {
  commandLine P(argc,argv," [-s] <inFile>");
  char* iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");
  bool compressed = P.getOptionValue("-c");
  bool binary = P.getOptionValue("-b");
  bool mmap = P.getOptionValue("-m");
  //cout << "mmap = " << mmap << endl;
  long rounds = P.getOptionLongValue("-rounds",1);

  size_t n_high_deg = P.getOptionLongValue("-nhighdeg", 4);

  string queryFileName = string(P.getOptionValue("-qf", ""));
  int combination_max = P.getOptionLongValue("-max_combination", 256);
  size_t bSize = P.getOptionLongValue("-batch", 4);

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
  std::shuffle(std::begin(userQueries), std::end(userQueries), rng);

  cout << "number of random queries: " << userQueries.size() << endl;
  int batch_size = userQueries.size();
  
  if (symmetric) {
    cout << "symmetric graph\n";
    graph<symmetricVertex> G =
      readGraph<symmetricVertex>(iFile,compressed,symmetric,binary,mmap); //symmetric graph
    // for(int r=0;r<rounds;r++) {
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

    unordered_set<intE> coverage;
    for (int i = 0; i < highdegQ.size(); i++) {
      intE* outnghs = G.V[highdegQ[i]].getOutNeighbors();
      long temp_degree =  G.V[highdegQ[i]].getOutDegree();
      for (long j = 0; j < temp_degree; j++) {
        coverage.insert(outnghs[j]);
      }
      cout << " top " << i+1 << " high degree has " << coverage.size() << " neighbors. " << 1.0*coverage.size() / n << endl;
    }

    // cout << "User queries degree distribution: \n";
    // for (int i = 0; i < userQueries.size(); i++) {
    //   long temp_degree =  G.V[userQueries[i]].getOutDegree();
    //   cout << temp_degree << endl;
    // }

    cout << "=================\n";
      
    // }
    G.del();
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

    unordered_set<intE> coverage;
    for (int i = 0; i < highdegQ.size(); i++) {
      long temp_degree =  G.V[highdegQ[i]].getOutDegree();
      cout << highdegQ[i] << ": " << temp_degree << endl;
      for (long j = 0; j < temp_degree; j++) {
        // cout << G.V[highdegQ[i]].getOutNeighbor(j) << endl;
        coverage.insert(G.V[highdegQ[i]].getOutNeighbor(j));
      }
      cout << " top " << i+1 << " high degree has " << coverage.size() << " neighbors. " << 1.0*coverage.size() / n << endl;
    }

    // cout << "User queries degree distribution: \n";
    // for (int i = 0; i < userQueries.size(); i++) {
    //   long temp_degree =  G.V[userQueries[i]].getOutDegree();
    //   cout << temp_degree << endl;
    // }
    cout << "=================\n";
      
    // }
    G.del();
  }
}

void test_6(int argc, char* argv[]) {
  commandLine P(argc,argv," [-s] <inFile>");
  char* iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");
  bool compressed = P.getOptionValue("-c");
  bool binary = P.getOptionValue("-b");
  bool mmap = P.getOptionValue("-m");
  //cout << "mmap = " << mmap << endl;
  long rounds = P.getOptionLongValue("-rounds",1);

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
  std::shuffle(std::begin(userQueries), std::end(userQueries), rng);

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
    for (int i = 0; i < combination_max; i++) {
      timer t_seq, t_batch, t_delay;
      std::shuffle(std::begin(batchedQuery), std::end(batchedQuery), rng);
      vector<long> tmp_batch;
      cout << "Evaluating queries: ";
      for (int j = 0; j < bSize; j++) {
        long tmp_query_id = batchedQuery[j];
        tmp_batch.push_back(tmp_query_id);
        cout << tmp_query_id << " ";
      }
      cout << endl;

      // Batching
      
      // Delayed batching
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
        cout << "No. " << j << " defer " << dist_to_high[j] << " iterations\n";
      }
      cout << "Total delays (delta): " << total_delays << endl;


      if (distances[tmp_batch[0]] == distances[tmp_batch[1]]) {
        cout << "to highest deg vtx: " << distances[tmp_batch[0]] << endl;
        Compute_Base(G,tmp_batch,P,true);

        uintE* query_distance;
        query_distance = Compute_Eval(G,tmp_batch,P);  // to get hops
        cout << "q0 to q1: " << query_distance[0+tmp_batch[1]*2] << endl;
        cout << "q1 to q0: " << query_distance[1+tmp_batch[0]*2] << endl;

        pbbs::delete_array(query_distance, n*tmp_batch.size());
      }
      

      // // To do: remove
      // for (int j = 0; j < tmp_batch.size(); j++) {
      //   vector<long> tmp_single_query;
      //   tmp_single_query.push_back(tmp_batch[j]);
      //   Compute_Base(G,tmp_single_query,P, true);
      // }
      cout << "=================\n";
    }
      
    // }
    G.del();
    pbbs::delete_array(distances, n);
    pbbs::delete_array(distances_multiple, n*highdegQ.size());
  }
}

// for dealying
void scenario3_fixed(int argc, char* argv[]) {
  commandLine P(argc,argv," [-s] <inFile>");
  char* iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");
  bool compressed = P.getOptionValue("-c");
  bool binary = P.getOptionValue("-b");
  bool mmap = P.getOptionValue("-m");
  //cout << "mmap = " << mmap << endl;
  long rounds = P.getOptionLongValue("-rounds",1);

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

  // // randomly shuffled each run
  // std::random_device rd;
  // auto rng = std::default_random_engine { rd() };
  // // auto rng = std::default_random_engine {};
  // std::shuffle(std::begin(userQueries), std::end(userQueries), rng);

  cout << "number of random queries: " << userQueries.size() << endl;
  int batch_size = userQueries.size();
  
  if (symmetric) {
    cout << "symmetric graph\n";
    graph<symmetricVertex> G =
      readGraph<symmetricVertex>(iFile,compressed,symmetric,binary,mmap); //symmetric graph
    // for(int r=0;r<rounds;r++) {
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
    for (int i = 0; i < high_deg_batch; i++) {
      highdegQ.push_back(vIDDegreePairs[i].first);
    }

    uintE* distances_multiple;
    distances_multiple = Compute_Eval(G,highdegQ,P);
    uintE* distances = pbbs::new_array<uintE>(n);
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
    int query_length = batchedQuery.size();
    vector<int> batch_size_array; // = {2, 4, 8, 16, 32, 64};
    batch_size_array.push_back(bSize);
    vector<double> seq_eval_times;
    vector<double> batch_eval_time;
    vector<double> delay_eval_time;
    // vector<double> b2_eval_times;
    // vector<double> b4_eval_times;
    // vector<double> b8_eval_times;
    // vector<double> b16_eval_times;
    // vector<double> b32_eval_times;
    // vector<double> b64_eval_times;
    cout << "=================\n";
    int start_from = P.getOptionLongValue("-start", 0);

    int tmp_batch_size = bSize;
    if (bSize == 1) {
        for (int i = start_from; i < query_length; i++) {
          timer t_seq; t_seq.start();
          vector<long> tmp_single_query;
          tmp_single_query.push_back(batchedQuery[i]);
          Compute_Base(G,tmp_single_query,P);
          t_seq.stop();
          double seq_time = t_seq.totalTime;
          cout << batchedQuery[i] << ": " << seq_time << endl;
          seq_eval_times.push_back(seq_time);
        }
    } else {
      int pos = start_from*bSize;
      for (; pos < query_length; pos+=bSize) {
        timer t_seq, t_batch, t_delay;
        vector<long> tmp_batch;
        cout << "Evaluating queries: ";
        for (int j = 0; j < bSize; j++) {
          long tmp_query_id = batchedQuery[pos+j];
          tmp_batch.push_back(tmp_query_id);
          cout << tmp_query_id << " ";
        }
        cout << endl;
        // Batching
        t_batch.start();
        Compute_Base(G,tmp_batch,P);
        t_batch.stop();

        // Delayed batching
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
          cout << "No. " << j << " defer " << dist_to_high[j] << " iterations\n";
        }
        cout << "Total delays (delta): " << total_delays << endl;
        t_delay.start();
        Compute_Delay(G,tmp_batch,P,dist_to_high);
        t_delay.stop();

        double batch_time = t_batch.totalTime;
        double delay_time = t_delay.totalTime;
        t_seq.reportTotal("sequential time");
        t_batch.reportTotal("batching evaluation time");
        t_delay.reportTotal("delayed batching evaluation time");

        batch_eval_time.push_back(batch_time);
        delay_eval_time.push_back(delay_time);
        // profiling the affinities
        Compute_Base(G,tmp_batch,P,true);
        Compute_Delay(G,tmp_batch,P,dist_to_high,true);

        cout << "=================\n";
      }
    }
    if (!seq_eval_times.empty()) {
      cout << "seq eval time: \n";
      for (int i = 0; i < seq_eval_times.size(); i++) {
        cout << seq_eval_times[i] << endl;
      }

      int k = 2;
      cout << "seq eval group by 2: \n";
      for (int i = 0; i < seq_eval_times.size(); i=i+k) {
        double tmp_sum = 0.0;
        for (int j = 0; j < k; j++) {
          tmp_sum += seq_eval_times[i+j];
        }
        cout << tmp_sum << endl;
      }

      k = 4;
      cout << "seq eval group by 4: \n";
      for (int i = 0; i < seq_eval_times.size(); i=i+k) {
        double tmp_sum = 0.0;
        for (int j = 0; j < k; j++) {
          tmp_sum += seq_eval_times[i+j];
        }
        cout << tmp_sum << endl;
      }

      k = 8;
      cout << "seq eval group by 8: \n";
      for (int i = 0; i < seq_eval_times.size(); i=i+k) {
        double tmp_sum = 0.0;
        for (int j = 0; j < k; j++) {
          tmp_sum += seq_eval_times[i+j];
        }
        cout << tmp_sum << endl;
      }

      k = 16;
      cout << "seq eval group by 16: \n";
      for (int i = 0; i < seq_eval_times.size(); i=i+k) {
        double tmp_sum = 0.0;
        for (int j = 0; j < k; j++) {
          tmp_sum += seq_eval_times[i+j];
        }
        cout << tmp_sum << endl;
      }

      k = 32;
      cout << "seq eval group by 32: \n";
      for (int i = 0; i < seq_eval_times.size(); i=i+k) {
        double tmp_sum = 0.0;
        for (int j = 0; j < k; j++) {
          tmp_sum += seq_eval_times[i+j];
        }
        cout << tmp_sum << endl;
      }

      k = 64;
      cout << "seq eval group by 64: \n";
      for (int i = 0; i < seq_eval_times.size(); i=i+k) {
        double tmp_sum = 0.0;
        for (int j = 0; j < k; j++) {
          tmp_sum += seq_eval_times[i+j];
        }
        cout << tmp_sum << endl;
      }
    }

    if (!batch_eval_time.empty()) {
      cout << "batch eval time: \n";
      for (int i = 0; i < batch_eval_time.size(); i++) {
        cout << batch_eval_time[i] << endl;
      }
    }

    if (!delay_eval_time.empty()) {
      cout << "delay eval time: \n";
      for (int i = 0; i < delay_eval_time.size(); i++) {
        cout << delay_eval_time[i] << endl;
      }
    }

    G.del();
    pbbs::delete_array(distances, n);
    pbbs::delete_array(distances_multiple, n*highdegQ.size());

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
    for (int i = 0; i < high_deg_batch; i++) {
      highdegQ.push_back(vIDDegreePairs[i].first);
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
    int query_length = batchedQuery.size();
    vector<int> batch_size_array; // = {2, 4, 8, 16, 32, 64};
    batch_size_array.push_back(bSize);
    vector<double> seq_eval_times;
    vector<double> batch_eval_time;
    vector<double> delay_eval_time;
    cout << "=================\n";
    int start_from = P.getOptionLongValue("-start", 0);
    cout << "start from: " << start_from << endl;
    int tmp_batch_size = bSize;
    if (bSize == 1) {
      cout << "Sequential eval query: \n";
        for (int i = start_from; i < query_length; i++) {
          timer t_seq; t_seq.start();
          vector<long> tmp_single_query;
          tmp_single_query.push_back(batchedQuery[i]);
          Compute_Base(G,tmp_single_query,P);
          t_seq.stop();
          double seq_time = t_seq.totalTime;
          cout << batchedQuery[i] << ": " << seq_time << endl;
          seq_eval_times.push_back(seq_time);
        }
    } else {
      int pos = start_from*bSize;
      for (; pos < query_length; pos+=bSize) {
        timer t_seq, t_batch, t_delay;
        vector<long> tmp_batch;
        cout << "Evaluating queries: ";
        for (int j = 0; j < bSize; j++) {
          long tmp_query_id = batchedQuery[pos+j];
          tmp_batch.push_back(tmp_query_id);
          cout << tmp_query_id << " ";
        }
        cout << endl;
        // Batching
        t_batch.start();
        Compute_Base(G,tmp_batch,P);
        t_batch.stop();

        // Delayed batching
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
          cout << "No. " << j << " defer " << dist_to_high[j] << " iterations\n";
        }
        cout << "Total delays (delta): " << total_delays << endl;
        t_delay.start();
        Compute_Delay(G,tmp_batch,P,dist_to_high);
        t_delay.stop();

        double batch_time = t_batch.totalTime;
        double delay_time = t_delay.totalTime;
        t_seq.reportTotal("sequential time");
        t_batch.reportTotal("batching evaluation time");
        t_delay.reportTotal("delayed batching evaluation time");

        batch_eval_time.push_back(batch_time);
        delay_eval_time.push_back(delay_time);
        // profiling the affinities
        Compute_Base(G,tmp_batch,P,true);
        Compute_Delay(G,tmp_batch,P,dist_to_high,true);
        cout << "=================\n";
      }
    }
    if (!seq_eval_times.empty()) {
      cout << "seq eval time: \n";
      for (int i = 0; i < seq_eval_times.size(); i++) {
        cout << seq_eval_times[i] << endl;
      }
      int k = 2;
      cout << "seq eval group by 2: \n";
      for (int i = 0; i < seq_eval_times.size(); i=i+k) {
        double tmp_sum = 0.0;
        for (int j = 0; j < k; j++) {
          tmp_sum += seq_eval_times[i+j];
        }
        cout << tmp_sum << endl;
      }

      k = 4;
      cout << "seq eval group by 4: \n";
      for (int i = 0; i < seq_eval_times.size(); i=i+k) {
        double tmp_sum = 0.0;
        for (int j = 0; j < k; j++) {
          tmp_sum += seq_eval_times[i+j];
        }
        cout << tmp_sum << endl;
      }

      k = 8;
      cout << "seq eval group by 8: \n";
      for (int i = 0; i < seq_eval_times.size(); i=i+k) {
        double tmp_sum = 0.0;
        for (int j = 0; j < k; j++) {
          tmp_sum += seq_eval_times[i+j];
        }
        cout << tmp_sum << endl;
      }

      k = 16;
      cout << "seq eval group by 16: \n";
      for (int i = 0; i < seq_eval_times.size(); i=i+k) {
        double tmp_sum = 0.0;
        for (int j = 0; j < k; j++) {
          tmp_sum += seq_eval_times[i+j];
        }
        cout << tmp_sum << endl;
      }

      k = 32;
      cout << "seq eval group by 32: \n";
      for (int i = 0; i < seq_eval_times.size(); i=i+k) {
        double tmp_sum = 0.0;
        for (int j = 0; j < k; j++) {
          tmp_sum += seq_eval_times[i+j];
        }
        cout << tmp_sum << endl;
      }

      k = 64;
      cout << "seq eval group by 64: \n";
      for (int i = 0; i < seq_eval_times.size(); i=i+k) {
        double tmp_sum = 0.0;
        for (int j = 0; j < k; j++) {
          tmp_sum += seq_eval_times[i+j];
        }
        cout << tmp_sum << endl;
      }
    }

    if (!batch_eval_time.empty()) {
      cout << "batch eval time: \n";
      for (int i = 0; i < batch_eval_time.size(); i++) {
        cout << batch_eval_time[i] << endl;
      }
    }

    if (!delay_eval_time.empty()) {
      cout << "delay eval time: \n";
      for (int i = 0; i < delay_eval_time.size(); i++) {
        cout << delay_eval_time[i] << endl;
      }
    }
      
    // }
    G.del();
    pbbs::delete_array(distances, n);
    pbbs::delete_array(distances_multiple, n*highdegQ.size());
  }
}

template <class vertex>
vector<long> reorderingByProperty(graph<vertex>& G, vector<long> truncatedQueries, int n_high_deg, commandLine P) {
  vector<long> sortedQueries;

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
  uintE* distances_multiple;
  G.transpose();
  cout << "sorting based on property values...\n";
  distances_multiple = Compute_Eval_Prop(G,highdegQ,P);
  uintE* distances = pbbs::new_array<uintE>(n);
  G.transpose();
  parallel_for(size_t i = 0; i < n; i++) {
    distances[i] = (uintE)MAXLEVEL;
  }
  parallel_for(size_t i = 0; i < n; i++) {
    for (int j = 0; j < high_deg_batch; j++) {
      if (distances_multiple[j+i*high_deg_batch] < distances[i]) {  // Note: only true for BFS and SSSP property.
        distances[i] = distances_multiple[j+i*high_deg_batch];
      }
    }
  }

  std::vector<std::pair<size_t, long>> vtxValuePairs;
  for (long i = 0; i < truncatedQueries.size(); i++) {
    vtxValuePairs.push_back(std::make_pair(truncatedQueries[i], distances[truncatedQueries[i]]));
  }

  std::sort(vtxValuePairs.begin(), vtxValuePairs.end(), sortByLargerSecondElement);
  for (int i = 0; i < vtxValuePairs.size(); i++) {
    long vID = vtxValuePairs[i].first;
    sortedQueries.push_back(vID);
    // cout << vID << " dist: " << degBatchPairs[i].second << endl;
  }

  cout << "\nsorted queries based on property values: \n";
  for (long i = 0; i < sortedQueries.size(); i++) {
    cout << sortedQueries[i] << ": " << distances[sortedQueries[i]] << endl;
  }

  return sortedQueries;
}

template <class vertex>
pair<vector<long>, vector<long>> streamingPreprocessing(graph<vertex>& G, vector<long> userQueries, int n_high_deg, int combination_max, commandLine P) {
  // input: graph, queries, queries_to_process, number of high degree vtxs.
  // return: unsorted and sorted queries of size queries_to_process
  vector<long> truncatedQueries;
  vector<long> sortedQueries;
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

  uintE* distances_multiple;
  // On edge reversed graph...
  G.transpose();
  distances_multiple = Compute_Eval(G,highdegQ,P);  // to get hops (bfs distance)
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

  for (long i = 0; i < combination_max; i++) {
    truncatedQueries.push_back(userQueries[i]);
  }
  std::vector<std::pair<size_t, long>> vtxHopPairs;

  cout << "\ninput queries: \n";
  for (long i = 0; i < truncatedQueries.size(); i++) {
    cout << truncatedQueries[i] << ": " << distances[truncatedQueries[i]] << endl;
    vtxHopPairs.push_back(std::make_pair(truncatedQueries[i], distances[truncatedQueries[i]]));
  }


  // int window_size = 256;
  // for (int i = 0; i < degBatchPairs.size()/window_size; i++) {
  //   std::sort(degBatchPairs.begin()+(i*window_size), degBatchPairs.begin()+(i+1)*window_size, sortByLargerSecondElement);
  // }

  std::sort(vtxHopPairs.begin(), vtxHopPairs.end(), sortByLargerSecondElement);

  for (int i = 0; i < vtxHopPairs.size(); i++) {
    long vID = vtxHopPairs[i].first;
    sortedQueries.push_back(vID);
    // cout << vID << " dist: " << degBatchPairs[i].second << endl;
  }

  cout << "\nsorted queries: \n";
  for (long i = 0; i < sortedQueries.size(); i++) {
    cout << sortedQueries[i] << ": " << distances[sortedQueries[i]] << endl;
  }

  return make_pair(truncatedQueries, sortedQueries);
}

template <class vertex>
void bufferStreaming(graph<vertex>& G, std::vector<long> bufferedQueries, int bSize, commandLine P, bool should_profile=false) {
  vector<double> latency_map;
    std::vector<double> arrivalTimes;
    for (int i = 0; i < bufferedQueries.size(); i++) {
      arrivalTimes.push_back(0.0);  // assuming all queries arrived at the same time.
    } 
    double static_latency = 0.0;
    // double start_time1 = update_timer.get_time();
    double earlier_start1 = 0;
    double new_est = 0;
    // double new_est = arrivalTimes[0+bSize-1];
    timer start_time1; start_time1.start();
    for (int i = 0; i < bufferedQueries.size(); i=i+bSize) {
      cout << "i: " << i << endl;
      std::vector<long> tmpBatch;
      double arrival_last_in_the_batch = arrivalTimes[i+bSize-1]; // last arrived in the batch
      for (int j = 0; j < bSize; j++) {
        tmpBatch.push_back(bufferedQueries[i+j]);
        if (arrival_last_in_the_batch < new_est) {
          // cout << "new_est - arrivalTimes[i+j]: " << new_est - arrivalTimes[i+j] << endl;
          static_latency += new_est - arrivalTimes[i+j];
        } else {
          static_latency += arrival_last_in_the_batch - arrivalTimes[i+j];
        }
      }

      timer t_t1;
      t_t1.start();
      Compute_Base(G,tmpBatch,P);
      t_t1.stop();
      double time1 = t_t1.totalTime;

      // record latency for each query
      for (int ii = 0; ii < bSize; ii++) {
        latency_map.push_back(new_est+time1);
      }

      if (arrival_last_in_the_batch < new_est) {
        new_est = time1 + new_est;
      } else {
        new_est = time1 + new_est + arrival_last_in_the_batch - new_est;
      }
      static_latency += (time1)*bSize;
      
      cout << "current latency: " << static_latency << endl;
    }
    start_time1.stop();
    double query_time1 = start_time1.totalTime;
    cout << "static batching version query time: " << query_time1 << endl;
    cout << "Static total latency: " << static_latency << endl;

    double check_sum = 0.0;
    sort(latency_map.begin(), latency_map.end());
    for (int ii = 0; ii < bufferedQueries.size(); ii++) {
      cout << latency_map[ii] << endl;
      check_sum += latency_map[ii];
    }
    cout << "check_sum: " << check_sum << endl;

    if (should_profile) {
      for (int i = 0; i < bufferedQueries.size(); i=i+bSize) {
        std::vector<long> tmpBatch;
        for (int j = 0; j < bSize; j++) {
          tmpBatch.push_back(bufferedQueries[i+j]);
        }
        Compute_Base(G,tmpBatch,P,true);
      }
    }
}

// for reordering
void scenario2(int argc, char* argv[]) {
  commandLine P(argc,argv," [-s] <inFile>");
  char* iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");
  bool compressed = P.getOptionValue("-c");
  bool binary = P.getOptionValue("-b");
  bool mmap = P.getOptionValue("-m");
  //cout << "mmap = " << mmap << endl;
  long rounds = P.getOptionLongValue("-rounds",1);

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
    long selection = P.getOptionLongValue("-order",1);
    if (selection == 1) {
      cout << "\nsequential evaluation..\n";
      bufferStreaming(G, truncatedQueries, 1, P);
    }
    if (selection == 2) {
      cout << "\non the unsorted buffer..\n";
      bufferStreaming(G, truncatedQueries, bSize, P, true);
    }
    if (selection == 3) {
      cout << "\non the sorted buffer..\n";
      bufferStreaming(G, sortedQueries, bSize, P, true);
    }
    if (selection == 4) {
      cout << "\non the property-based sorted buffer..\n";
      propertySortedQueries = reorderingByProperty(G, truncatedQueries, n_high_deg, P);
      bufferStreaming(G, propertySortedQueries, bSize, P, true);
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
    tie(truncatedQueries, sortedQueries) = streamingPreprocessing(G, userQueries, n_high_deg, combination_max, P);
    
    // start streaming.
    // input: G, P, bufferedQueries, batch size
    long selection = P.getOptionLongValue("-order",1);
    if (selection == 1) {
      cout << "\nsequential evaluation..\n";
      bufferStreaming(G, truncatedQueries, 1, P);
    }
    if (selection == 2) {
      cout << "\non the unsorted buffer..\n";
      bufferStreaming(G, truncatedQueries, bSize, P, true);
    }
    if (selection == 3) {
      cout << "\non the sorted buffer..\n";
      bufferStreaming(G, sortedQueries, bSize, P, true);
    }
    if (selection == 4) {
      cout << "\non the property-based sorted buffer..\n";
      propertySortedQueries = reorderingByProperty(G, truncatedQueries, n_high_deg, P);
      bufferStreaming(G, propertySortedQueries, bSize, P, true);
    }
    // cout << "\non the property-based sorted buffer..\n";
    // bufferStreaming(G, propertySortedQueries, bSize, P, true);

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
                // std::cout << " " << i+1;
            }
        }
        ret.push_back(tmp_sample);
        if (ret.size() >= max_size) {   // break early
            return ret;
        }
        // std::cout << std::endl;
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
    return ret;
}

// for reordering
void test_3(int argc, char* argv[]) {
  commandLine P(argc,argv," [-s] <inFile>");
  char* iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");
  bool compressed = P.getOptionValue("-c");
  bool binary = P.getOptionValue("-b");
  bool mmap = P.getOptionValue("-m");
  //cout << "mmap = " << mmap << endl;
  long rounds = P.getOptionLongValue("-rounds",1);

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
  std::shuffle(std::begin(userQueries), std::end(userQueries), rng);
  cout << "number of random queries: " << userQueries.size() << endl;

  int setSize = userQueries.size();
  std::vector<long> testQuery[setSize];

  if (symmetric) {
    cout << "symmetric graph: not implemented!\n";

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

    vector<long> sorted_first;
    vector<long> sorted_second;
    vector<long> unsorted_first;
    vector<long> unsorted_second;

    int buffer_size = 16;
    vector<long> tmp_queries_first;
    vector<long> tmp_queries_second;
    for (int i = 0; i < buffer_size/2; i++) {
      tmp_queries_first.push_back(userQueries[i]);
      tmp_queries_second.push_back(userQueries[i+buffer_size/2]);
    }

    // tie(unsorted_first, sorted_first) = streamingPreprocessing(G, tmp_queries_first, n_high_deg, buffer_size/2, P);
    // tie(unsorted_second, sorted_second) = streamingPreprocessing(G, tmp_queries_second, n_high_deg, buffer_size/2, P);

    // for (int i = 0; i < buffer_size/2; i++) {
    //   sortedQueries.push_back(sorted_first[i]);
    //   truncatedQueries.push_back(unsorted_first[i]);
    // }
    // for (int i = 0; i < buffer_size/2; i++) {
    //   sortedQueries.push_back(sorted_second[i]);
    //   truncatedQueries.push_back(unsorted_second[i]);
    // }

    tie(truncatedQueries, sortedQueries) = streamingPreprocessing(G, userQueries, n_high_deg, buffer_size, P);

    // start streaming.
    // input: G, P, bufferedQueries, batch size
    long selection = P.getOptionLongValue("-order",1);
    if (selection == 1) {
      cout << "\nsequential evaluation..\n";
      bufferStreaming(G, truncatedQueries, 1, P);
    }
    if (selection == 2) {
      cout << "\non the unsorted buffer..\n";
      bufferStreaming(G, truncatedQueries, bSize, P, false);
    }
    if (selection == 3) {
      cout << "\non the sorted buffer..\n";
      bufferStreaming(G, sortedQueries, bSize, P, false);
    }
    if (selection == 4) {
      cout << "\nenumeration..\n";
      
      vector<vector<long>> combs = my_comb(16, 8);
      for (int j = 0; j < combs.size(); j++) {
        vector<long> tmp_queries;
        vector<int> remains;
        cout << j << ": ";
        for (int k = 0; k < combs[j].size(); k++) {
          cout << combs[j][k] << " ";
          remains.push_back(userQueries.size()-1-combs[j][k]);
          tmp_queries.push_back(userQueries[combs[j][k]]);
        }
        for (int k = 0; k < remains.size(); k++) {
          cout << remains[k] << " ";
          tmp_queries.push_back(userQueries[remains[k]]);
        }
        cout << endl;
        bufferStreaming(G, tmp_queries, bSize, P, false);
      }
      // bufferStreaming(G, sortedQueries, bSize, P, false);
    }
    // cout << "\non the property-based sorted buffer..\n";
    // bufferStreaming(G, propertySortedQueries, bSize, P, true);

  }
}

void test_1(int argc, char* argv[], vector<vector<long>> my_comb) {
  commandLine P(argc,argv," [-s] <inFile>");
  char* iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");
  bool compressed = P.getOptionValue("-c");
  bool binary = P.getOptionValue("-b");
  bool mmap = P.getOptionValue("-m");
  //cout << "mmap = " << mmap << endl;
  long rounds = P.getOptionLongValue("-rounds",1);

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
  std::shuffle(std::begin(userQueries), std::end(userQueries), rng);

  cout << "number of random queries: " << userQueries.size() << endl;
  
  if (symmetric) {
    cout << "symmetric graph\n";
    graph<symmetricVertex> G =
      readGraph<symmetricVertex>(iFile,compressed,symmetric,binary,mmap); //symmetric graph
    // for(int r=0;r<rounds;r++) {
    cout << "n=" << G.n << " m=" << G.m << endl;
    size_t n = G.n;

    // average (random) batch evaluation
    for (int i = 0; i < userQueries.size(); i+=bSize) {
      timer t_seq, t_batch, t_delay;
      vector<long> tmp_batch;
      for (int j = 0; j < bSize; j++) {
        long tmp_query_id = userQueries[i+j];
        tmp_batch.push_back(tmp_query_id);
        cout << tmp_query_id << " ";
      }
      cout << endl;

      // Sequential
      t_seq.start();
      for (int j = 0; j < tmp_batch.size(); j++) {
        vector<long> tmp_single_query;
        tmp_single_query.push_back(tmp_batch[j]);
        Compute_Base(G,tmp_single_query,P);
      }
      t_seq.stop();

      // Batching
      t_batch.start();
      // Compute_Base(G,tmp_batch,P);
      t_batch.stop();

      double seq_time = t_seq.totalTime;
      double batch_time = t_batch.totalTime / rounds;
      t_seq.reportTotal("random sequential time");
      // t_batch.reportTotal("random batching evaluation time");

      // cout << "random batching speedup: " << seq_time / batch_time << endl;

      // // profiling the affinities
      // Compute_Base(G,tmp_batch,P,true);

      cout << "=================\n";
      
    }
    
    cout << "Exhausive evaluation...\n";
    
    // Query evaluation: sequential, batching, delayed batching
    vector<long> batchedQuery;
    batchedQuery = userQueries;
    cout << "=================\n";
    for (int i = 0; i < my_comb.size(); i++) {
      timer t_seq, t_batch, t_delay;
      vector<long> tmp_batch;
      cout << "Evaluating queries: ";
      auto tmp_comb = my_comb[i];
      for (int j = 0; j < bSize; j++) {
        long tmp_query_id = batchedQuery[tmp_comb[j]];
        tmp_batch.push_back(tmp_query_id);
        cout << tmp_query_id << " ";
      }
      cout << endl;

      // // Sequential
      // t_seq.start();
      // for (int j = 0; j < tmp_batch.size(); j++) {
      //   vector<long> tmp_single_query;
      //   tmp_single_query.push_back(tmp_batch[j]);
      //   Compute_Base(G,tmp_single_query,P);
      // }
      // t_seq.stop();

      // Batching
      t_batch.start();
      Compute_Base(G,tmp_batch,P);
      t_batch.stop();

      // double seq_time = t_seq.totalTime;
      double batch_time = t_batch.totalTime / rounds;
      // t_seq.reportTotal("sequential time");
      t_batch.reportTotal("batching evaluation time");

      // cout << "Batching speedup: " << seq_time / batch_time << endl;

      // // profiling the affinities
      // Compute_Base(G,tmp_batch,P,true);

      cout << "=================\n";
    }
    G.del();

  } else {
    // For directed graph...
    cout << "asymmetric graph\n";
    graph<asymmetricVertex> G =
      readGraph<asymmetricVertex>(iFile,compressed,symmetric,binary,mmap); //asymmetric graph
    cout << "n=" << G.n << " m=" << G.m << endl;

    size_t n = G.n;
    // ========================================
    
    // average (random) batch evaluation
    for (int i = 0; i < userQueries.size(); i+=bSize) {
      timer t_seq, t_batch, t_delay;
      vector<long> tmp_batch;
      for (int j = 0; j < bSize; j++) {
        long tmp_query_id = userQueries[i+j];
        tmp_batch.push_back(tmp_query_id);
        cout << tmp_query_id << " ";
      }
      cout << endl;

      // Sequential
      t_seq.start();
      for (int j = 0; j < tmp_batch.size(); j++) {
        vector<long> tmp_single_query;
        tmp_single_query.push_back(tmp_batch[j]);
        Compute_Base(G,tmp_single_query,P);
      }
      t_seq.stop();

      // Batching
      t_batch.start();
      // Compute_Base(G,tmp_batch,P);
      t_batch.stop();

      double seq_time = t_seq.totalTime;
      double batch_time = t_batch.totalTime / rounds;
      t_seq.reportTotal("random sequential time");
      // t_batch.reportTotal("random batching evaluation time");

      // cout << "random batching speedup: " << seq_time / batch_time << endl;

      // // profiling the affinities
      // Compute_Base(G,tmp_batch,P,true);

      cout << "=================\n";
      
    }
    
    cout << "Exhausive evaluation...\n";
    
    // Query evaluation: sequential, batching, delayed batching
    vector<long> batchedQuery;
    batchedQuery = userQueries;
    cout << "=================\n";
    for (int i = 0; i < my_comb.size(); i++) {
      timer t_seq, t_batch, t_delay;
      vector<long> tmp_batch;
      cout << "Evaluating queries: ";
      auto tmp_comb = my_comb[i];
      for (int j = 0; j < bSize; j++) {
        long tmp_query_id = batchedQuery[tmp_comb[j]];
        tmp_batch.push_back(tmp_query_id);
        cout << tmp_query_id << " ";
      }
      cout << endl;

      // // Sequential
      // t_seq.start();
      // for (int j = 0; j < tmp_batch.size(); j++) {
      //   vector<long> tmp_single_query;
      //   tmp_single_query.push_back(tmp_batch[j]);
      //   Compute_Base(G,tmp_single_query,P);
      // }
      // t_seq.stop();

      // Batching
      t_batch.start();
      Compute_Base(G,tmp_batch,P);
      t_batch.stop();

      // double seq_time = t_seq.totalTime;
      double batch_time = t_batch.totalTime / rounds;
      // t_seq.reportTotal("sequential time");
      t_batch.reportTotal("batching evaluation time");

      // cout << "Batching speedup: " << seq_time / batch_time << endl;

      // // profiling the affinities
      // Compute_Base(G,tmp_batch,P,true);

      cout << "=================\n";
    }
      
    // }
    G.del();
  }
}

// for dealying
void test_2(int argc, char* argv[]) {
  commandLine P(argc,argv," [-s] <inFile>");
  char* iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");
  bool compressed = P.getOptionValue("-c");
  bool binary = P.getOptionValue("-b");
  bool mmap = P.getOptionValue("-m");
  //cout << "mmap = " << mmap << endl;
  long rounds = P.getOptionLongValue("-rounds",1);

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
  std::shuffle(std::begin(userQueries), std::end(userQueries), rng);

  cout << "number of random queries: " << userQueries.size() << endl;
  int batch_size = userQueries.size();
  
  if (symmetric) {
    cout << "symmetric graph not implemented\n";

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
    for (int i = 0; i < high_deg_batch; i++) {
      highdegQ.push_back(vIDDegreePairs[i].first);
    }
    uintE* distances_multiple;
    G.transpose();
    cout << "sorting based on property values...\n";
    distances_multiple = Compute_Eval_Prop(G,highdegQ,P);
    uintE* distances = pbbs::new_array<uintE>(n);
    G.transpose();
    parallel_for(size_t i = 0; i < n; i++) {
      distances[i] = (uintE)MAXLEVEL;
    }
    parallel_for(size_t i = 0; i < n; i++) {
      for (int j = 0; j < high_deg_batch; j++) {
        if (distances_multiple[j+i*high_deg_batch] < distances[i]) {  // Note: only true for BFS and SSSP property.
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
    for (int i = 0; i < combination_max; i++) {
      timer t_seq, t_batch, t_delay;
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
      t_seq.start();
      for (int j = 0; j < tmp_batch.size(); j++) {
        vector<long> tmp_single_query;
        tmp_single_query.push_back(tmp_batch[j]);
        Compute_Base(G,tmp_single_query,P);
      }
      t_seq.stop();

      // Batching
      t_batch.start();
      Compute_Base(G,tmp_batch,P);
      t_batch.stop();

      // Delayed batching
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
        cout << "No. " << j << " defer " << dist_to_high[j] << " iterations\n";
      }
      cout << "Total delays (delta): " << total_delays << endl;


      double seq_time = t_seq.totalTime;
      double batch_time = t_batch.totalTime / rounds;
      double delay_time = t_delay.totalTime / rounds;
      t_seq.reportTotal("sequential time");
      t_batch.reportTotal("batching evaluation time");
      t_delay.reportTotal("delayed batching evaluation time");

      cout << "Batching speedup: " << seq_time / batch_time << endl;
      cout << "Delayed batching speedup: " << seq_time / delay_time << endl;

      // profiling the affinities
      Compute_Base(G,tmp_batch,P,true);

    }
      
    // }
    G.del();
    pbbs::delete_array(distances, n);
    pbbs::delete_array(distances_multiple, n*highdegQ.size());
  }
}

int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv," [-s] <inFile>");
  string options = string(P.getOptionValue("-option", "scenario3"));
  if (options == "scenario3") {
    cout << "running scenraio 3\n";
    scenario3(argc, argv); // delaying
  }
  if (options == "scenario3_fixed") {
    cout << "running scenraio 3\n";
    scenario3_fixed(argc, argv); // delaying
  }
    
  if (options == "scenario2") {
    cout << "running scenraio 2\n";
    scenario2(argc, argv);  // reordering
  }

  if (options == "adaptive") {
    cout << "adaptive batching...\n";
    scenario_adaptive(argc, argv);
  }

  if (options == "adaptive_new") {
    cout << "adaptive batching (new)...\n";
    scenario_adaptive_new(argc, argv);
  }
    
  if (options == "random_generation") {
    cout << "random query generation\n";
    QueryGeneration_Random(argc, argv);
  }

  if (options == "test_1") {
    cout << "testing 16 queries exhaustively\n";
    size_t bSize = P.getOptionLongValue("-batch", 4);
    auto test_comb = my_comb(16, (int)bSize);
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(test_comb.begin(), test_comb.end(), g);
    cout << "test combinations: " << test_comb.size() << endl;
    test_1(argc, argv, test_comb);
  }

  if (options == "test_2") {
    cout << "testing property values and affinity\n";
    test_2(argc, argv);
  }

  if (options == "test_3") {
    cout << "testing reordering with the best speedup\n";
    test_3(argc, argv);
  }

  if (options == "test_4") {
    cout << "testing delaying with the best speedup\n";
    test_4(argc, argv);
  }

  if (options == "test_5") {
    cout << "testing high degree vertices coverage\n";
    test_5(argc, argv);
  }

  if (options == "test_6") {
    cout << "testing hops between two queries and hops to the high degree\n";
    test_6(argc, argv);
  }
    
}
#endif