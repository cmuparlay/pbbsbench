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

#define NOTMAIN 1
#include <iostream>
#include <limits.h>
#include "parlay/primitives.h"
#include "parlay/parallel.h"
#include "common/graph.h"
#include "common/speculative_for.h"
#include "algorithm/union_find.h"
#include "ST.h"

using reservation = pbbs::reservation<edgeId>;

struct unionFindStep {
  vertexId u;  vertexId v;  
  edgeArray<vertexId> const &E;
  parlay::sequence<reservation> &R;
  unionFind<vertexId> &UF;
  unionFindStep(edgeArray<vertexId> const &E,
		unionFind<vertexId> &UF,
		parlay::sequence<reservation> &R)
    : E(E), R(R), UF(UF) {} 

  bool reserve(edgeId i) {
    u = UF.find(E[i].u);
    v = UF.find(E[i].v);
    if (u > v) std::swap(u,v);
    if (u != v) {
      R[v].reserve(i);
      return 1;
    } else return 0;
  }

  bool commit(edgeId i) {
    if (R[v].check(i)) { UF.link(v, u); return 1; }
    else return 0;
  }
};

parlay::sequence<edgeId> st(edgeArray<vertexId> const &G){
  size_t m = G.nonZeros;
  size_t n = G.numRows;
  unionFind<vertexId> UF(n);
  parlay::sequence<reservation> R(n);
  unionFindStep UFStep(G, UF, R);
  pbbs::speculative_for<edgeId>(UFStep, 0, m, 100);
  return parlay::internal::filter_map(R,
		  [&] (const reservation& a) -> bool {return a.reserved();},
		  [&] (const reservation& a) -> edgeId {return a.get();});
} 
