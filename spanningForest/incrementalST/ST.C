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
#include "sequence.h"
#include "parallel.h"
#include "graph.h"
#include "ST.h"
#include "speculative_for.h"
#include "union_find.h"

template <class intT>
struct unionFindStep {
  intT u;  intT v;  
  edgeArray<intT> const &E;
  pbbs::sequence<reservation> &R;
  unionFind &UF;
  unionFindStep(edgeArray<intT> const &E,
		unionFind &UF,
		pbbs::sequence<reservation> &R)
    : E(E), R(R), UF(UF) {} 

  bool reserve(intT i) {
    u = UF.find(E[i].u);
    v = UF.find(E[i].v);
    if (u > v) std::swap(u,v);
    if (u != v) {
      R[v].reserve(i);
      return 1;
    } else return 0;
  }

  bool commit(intT i) {
    if (R[v].check(i)) { UF.link(v, u); return 1; }
    else return 0;
  }
};

std::pair<intT*,intT> st(edgeArray<intT> const &G){
  intT m = G.nonZeros;
  intT n = G.numRows;
  unionFind UF(n);
  pbbs::sequence<reservation> R(n);
  unionFindStep<intT> UFStep(G, UF, R);
  speculative_for(UFStep, 0, m, 100);
  auto stIdx = pbbs::filter(R, [&] (reservation a) {
      return a.reserved();});
  size_t n_out = stIdx.size();
  cout << "Tree size = " << n_out << endl;
  return std::pair<intT*,intT>((intT*) stIdx.to_array(), n_out);
}
