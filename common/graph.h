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
  edge(intT f, intT s) : u(f), v(s) {}
};

template <class intT>
struct edgeArray {
  edge<intT>* E;
  intT numRows;
  intT numCols;
  intT nonZeros;
  void del() {free(E);}
  edgeArray(edge<intT> *EE, intT r, intT c, intT nz) :
    E(EE), numRows(r), numCols(c), nonZeros(nz) {}
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
  wghEdge<intT> *E;
  intT n; intT m;
  wghEdgeArray(wghEdge<intT>* EE, intT nn, intT mm) : E(EE), n(nn), m(mm) {}
  void del() { free(E);}
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
struct graph {
  vertex<intT> *V;
  intT n;
  intT m;
  intT* allocatedInplace;
  graph(vertex<intT>* VV, intT nn, uintT mm) 
    : V(VV), n(nn), m(mm), allocatedInplace(NULL) {}
  graph(vertex<intT>* VV, intT nn, uintT mm, intT* ai) 
    : V(VV), n(nn), m(mm), allocatedInplace(ai) {}
  intT* vertices() { return allocatedInplace+2; }
  intT* edges() { return allocatedInplace+2+n; }
  graph copy() {
    vertex<intT>* VN = newA(vertex<intT>,n);
    intT* _allocatedInplace = newA(intT,n+m+2);
    _allocatedInplace[0] = n;
    _allocatedInplace[1] = m;
    intT* Edges = _allocatedInplace+n+2;
    intT k = 0;
    for (intT i=0; i < n; i++) {
      _allocatedInplace[i+2] = allocatedInplace[i+2];
      VN[i] = V[i];
      VN[i].Neighbors = Edges + k;
      for (intT j =0; j < V[i].degree; j++) 
	Edges[k++] = V[i].Neighbors[j];
    }
    return graph(VN, n, m, _allocatedInplace);
  } 
  void del() {
    if (allocatedInplace == NULL) 
      for (intT i=0; i < n; i++) V[i].del();
    else free(allocatedInplace);
    free(V);
  }
};

template <class intT, class intE>
struct graphC {
  long n,m;
  intT* offsets;
  intE* edges;
  graphC(intT* _offsets, intE* _edges, long _n, long _m)
  : offsets(_offsets), edges(_edges), n(_n), m(_m) {}
  graphC copy() {
    intT* o = newA(intT,n+1);
    intE* e = newA(intE,m);
    parallel_for(0, n+1, [&] (size_t i) {o[i] = offsets[i];});
    parallel_for(0, m, [&] (size_t i) {e[i] = edges[i];});
    return graphC(o,e,n,m);
  }
  size_t degree(size_t i) {
    return offsets[i+1] - offsets[i];
  }
  void del() {
    free(offsets); free(edges);
  }
};

template <class intT>
struct wghGraph {
  wghVertex<intT> *V;
  intT n;
  uintT m;
  intT* allocatedInplace;
  intT* weights;
wghGraph(wghVertex<intT>* VV, intT nn, uintT mm) 
    : V(VV), n(nn), m(mm), allocatedInplace(NULL) {}
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
    else { free(allocatedInplace); free(weights); }
    free(V);
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

