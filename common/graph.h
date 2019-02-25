// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
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

#pragma once

#include <iostream>
#include <algorithm>
#include "../pbbslib/parallel.h"
#include "glue.h"

// **************************************************************
//    SPARSE ROW MAJOR REPRESENTATION
// **************************************************************

template <class ETYPE, class intT>
struct sparseRowMajor {
  intT numRows;
  intT numCols;
  intT nonZeros;
  intT* Starts;
  intT* ColIds;
  ETYPE* Values;
  void del() {free(Starts); free(ColIds); if (Values != NULL) free(Values);}
  sparseRowMajor(intT n, intT m, intT nz, intT* S, intT* C, ETYPE* V) :
    numRows(n), numCols(m), nonZeros(nz), 
    Starts(S), ColIds(C), Values(V) {}
};

//typedef sparseRowMajor<double> sparseRowMajorD;

// **************************************************************
//    EDGE ARRAY REPRESENTATION
// **************************************************************

template <class intT>
struct edge {
  intT u;
  intT v;
  edge() {}
  edge(intT f, intT s) : u(f), v(s) {}
};

template <class intT>
struct edgeArray {
  pbbs::sequence<edge<intT>> E;
  size_t numRows;
  size_t numCols;
  size_t nonZeros;
  void del() {}
  edgeArray(pbbs::sequence<edge<intT>> EE, size_t r, size_t c) :
    E(std::move(EE)), numRows(r), numCols(c), nonZeros(E.size()) {}
  edgeArray() {}
};

// **************************************************************
//    WEIGHED EDGE ARRAY
// **************************************************************

template <class intT>
struct wghEdge {
  intT u, v;
  double weight;
  wghEdge() {}
  wghEdge(intT _u, intT _v, double w) : u(_u), v(_v), weight(w) {}
};

template <class intT>
struct wghEdgeArray {
  pbbs::sequence<wghEdge<intT>> E;
  intT n; intT m;
  wghEdgeArray(pbbs::sequence<wghEdge<intT>> E, intT n) : E(E), n(n) {}
  void del() {}
};

// **************************************************************
//    ADJACENCY ARRAY REPRESENTATION
// **************************************************************

template <class intT>
struct vertex {
  intT* Neighbors;
  intT degree;
  void del() {free(Neighbors);}
  vertex(intT* N, intT d) : Neighbors(N), degree(d) {}
  vertex() : Neighbors(NULL), degree(0) {}
};

template <class intT>
struct graph {
  pbbs::sequence<intT> offsets;
  pbbs::sequence<intT> edges;
  size_t n;
  size_t m;
  size_t numVertices() const {return n;}
  size_t numEdges() const {return m;}
  const pbbs::sequence<intT>& get_offsets() const {
    return offsets;
  }
  vertex<intT> operator[] (const size_t i) const {
    return vertex<intT>(edges.begin() + offsets[i],
			offsets[i+1] - offsets[i]);}
  
  graph(pbbs::sequence<intT> offsets_, pbbs::sequence<intT> edges_, size_t n) 
    : offsets(std::move(offsets_)), edges(std::move(edges_)), n(n), m(edges.size()) {
    if (offsets.size() != n + 1) { cout << "error in graph constructor" << endl;}
  }
};

template <class intT, class intE>
struct graphC {
  long n,m;
  pbbs::sequence<intT> offsets;
  pbbs::sequence<intE> edges;
  graphC(pbbs::sequence<intT> offsets, pbbs::sequence<intE> edges)
    : offsets(offsets), edges(edges), n(offsets.size()-1), m(edges.size()) {}
  graphC copy() {return graphC(offsets, edges);}
  size_t degree(size_t i) {return offsets[i+1] - offsets[i]; }
};

template <class intT>
struct wghVertex {
  intT* Neighbors;
  intT degree;
  intT* nghWeights;
  void del() {free(Neighbors); free(nghWeights);}
  wghVertex(intT* N, intT* W, intT d) : Neighbors(N), nghWeights(W), degree(d) {}
};

template <class intT>
struct wghGraph {
  wghVertex<intT> *V;
  intT n;
  uintT m;
  intT* allocatedInplace;
  intT* weights;
  wghGraph(wghVertex<intT>* VV, intT nn, uintT mm) 
    : V(VV), n(nn), m(mm), allocatedInplace(NULL) {
    cout << "double check correctness in wghGraph" << endl;}
  wghGraph(wghVertex<intT>* VV, intT nn, uintT mm, intT* ai, intT* _weights) 
    : V(VV), n(nn), m(mm), allocatedInplace(ai), weights(_weights) {}
  wghGraph copy() {
    wghVertex<intT>* VN = newA(wghVertex<intT>,n);
    intT* Edges = newA(intT,m);
    intT* Weights = newA(intT,m);
    intT k = 0;
    for (intT i=0; i < n; i++) {
      VN[i] = V[i];
      VN[i].Neighbors = Edges + k;
      VN[i].nghWeights = Weights + k;
      for (intT j =0; j < V[i].degree; j++){ 
	Edges[k] = V[i].Neighbors[j];
	Weights[k++] = V[i].nghWeights[j];
      }
    }
    return wghGraph(VN, n, m, Edges, Weights);
  } 
  void del() {
    if (allocatedInplace == NULL) 
      for (intT i=0; i < n; i++) V[i].del();
    else { pbbs::free_array(allocatedInplace); }
    pbbs::free_array(V);
  }
};

template <typename intT>
struct FlowGraph {
  wghGraph<intT> g;
  intT source, sink;
  FlowGraph(wghGraph<intT> g, intT source, intT sink)
    : g(g), source(source), sink(sink) {}
  FlowGraph copy() {
    return FlowGraph(g.copy(), source, sink);
  }
  void del() { g.del(); }
};

