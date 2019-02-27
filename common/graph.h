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
#include "../pbbslib/sequence.h"
using intT = int;
using uintT = unsigned int;

// **************************************************************
//    SPARSE ROW MAJOR REPRESENTATION
// **************************************************************

// template <class ETYPE, class intT>
// struct sparseRowMajor {
//   intT numRows;
//   intT numCols;
//   intT nonZeros;
//   intT* Starts;
//   intT* ColIds;
//   ETYPE* Values;
//   void del() {free(Starts); free(ColIds); if (Values != NULL) free(Values);}
//   sparseRowMajor(intT n, intT m, intT nz, intT* S, intT* C, ETYPE* V) :
//     numRows(n), numCols(m), nonZeros(nz), 
//     Starts(S), ColIds(C), Values(V) {}
// };

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
  edgeArray(pbbs::sequence<edge<intT>> EE, size_t r, size_t c) :
    E(std::move(EE)), numRows(r), numCols(c), nonZeros(E.size()) {}
  edgeArray() {}
  edge<intT> operator[] (const size_t i) const {return E[i];}
};

// **************************************************************
//    WEIGHED EDGE ARRAY
// **************************************************************

template <class intT>
struct wghEdge {
  intT u, v;
  double weight;
  wghEdge() {}
  wghEdge(intT _u, intT _v, float w) : u(_u), v(_v), weight(w) {}
};

template <class intT>
struct wghEdgeArray {
  pbbs::sequence<wghEdge<intT>> E;
  intT n; intT m;
  wghEdgeArray(pbbs::sequence<wghEdge<intT>> E_, intT n) 
    : E(std::move(E_)), n(n), m(E.size()) {}
  wghEdgeArray() {}
  wghEdge<intT> operator[] (const size_t i) const {return E[i];}
};

// **************************************************************
//    ADJACENCY ARRAY REPRESENTATION
// **************************************************************

template <class intT>
struct vertex {
  intT* Neighbors;
  intT degree;
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
  
  graph(pbbs::sequence<intT> offsets_,
	pbbs::sequence<intT> edges_,
	size_t n) 
    : offsets(std::move(offsets_)), edges(std::move(edges_)), n(n), m(edges.size()) {
    if (offsets.size() != n + 1) { cout << "error in graph constructor" << endl;}
  }
};

// **************************************************************
//    WEIGHTED ADJACENCY ARRAY REPRESENTATION
// **************************************************************

template <class intT>
struct wghVertex {
  intT* Neighbors;
  intT degree;
  intT* nghWeights;
  wghVertex(intT* N, intT* W, intT d) : Neighbors(N), nghWeights(W), degree(d) {}
};

template <class intT>
struct wghGraph {
  pbbs::sequence<intT> offsets;
  pbbs::sequence<intT> edges;
  pbbs::sequence<intT> weights;
  intT n;
  uintT m;
  size_t numVertices() const {return n;}
  size_t numEdges() const {return m;}
  const pbbs::sequence<intT>& get_offsets() const {
    return offsets;
  }
  wghVertex<intT> operator[] (const size_t i) const {
    return wghVertex<intT>(edges.begin() + offsets[i],
			   weights.begin() + offsets[i],
			   offsets[i+1] - offsets[i]);}
  wghGraph(pbbs::sequence<intT> offsets_,
	   pbbs::sequence<intT> edges_,
	   pbbs::sequence<intT> weights_,
	   size_t n) 
    : offsets(std::move(offsets_)), edges(std::move(edges_)),
      weights(std::move(weights_)), n(n), m(edges.size()) {
    if (offsets.size() != n + 1 || weights.size() != edges.size()) {
      cout << "error in weighted graph constructor" << endl;}
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

