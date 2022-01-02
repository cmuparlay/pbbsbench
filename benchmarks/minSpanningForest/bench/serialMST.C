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

#include <iostream>
#include <algorithm>
#include "parlay/parallel.h"
#include "common/graph.h"
#include "common/get_time.h"
#include "algorithm/union_find.h"
#include "MST.h"

using namespace std;

// **************************************************************
//    SERIAL MST USING KRUSKAL's ALGORITHM WITH FILTERING
// **************************************************************

struct edgeAndIndex {
  vertexId u;
  vertexId v;
  double weight;
  edgeId id;
  edgeAndIndex() {}
  edgeAndIndex(vertexId _u, vertexId _v, double w, edgeId _id) 
    : u(_u), v(_v), id(_id), weight(w) {}
};

int unionFindLoop(edgeAndIndex* E, size_t m, size_t nInMst,
		  unionFind<vertexId> &UF, parlay::sequence<edgeId> &mst) {
  for (size_t i = 0; i < m; i++) {
    vertexId u = UF.find(E[i].u);
    vertexId v = UF.find(E[i].v);
    
    // union by rank
    if (u != v) {
      UF.union_roots(u,v);
      mst.push_back(E[i].id); // add edge to output
    }
  }
  return nInMst;
}

parlay::sequence<edgeId> mst(wghEdgeArray<vertexId,edgeWeight> &G) {
  // tag with edge id
  auto EI = parlay::tabulate(G.m, [&] (size_t i) -> edgeAndIndex {
      return edgeAndIndex(G.E[i].u, G.E[i].v, G.E[i].weight, i);});

  auto edgeLess = [&] (edgeAndIndex a, edgeAndIndex b) {
    return (a.weight == b.weight) ? (a.id < b.id) : (a.weight < b.weight);};

  // partition so minimum 4/3 n elements are at bottom
  size_t l = min(4*G.n/3,G.m);
  std::nth_element(EI.begin(), EI.begin()+l, EI.begin()+G.m, edgeLess);

  // sort the prefix
  std::sort(EI.begin(), EI.begin()+l, edgeLess);

  // create union-find structure
  unionFind<vertexId> UF(G.n);

  // mst edges added to this sequence
  parlay::sequence<edgeId> mst;

  // run union-find over the prefix
  size_t nInMst = unionFindLoop(EI.begin(), l, 0, UF, mst);

  // pack down active edges
  size_t k = 0;
  for (size_t i = l; i < G.m; i++) {
    vertexId u = UF.find(EI[i].u);
    vertexId v = UF.find(EI[i].v);
    if (u != v) EI[l + k++] = EI[i]; 
  }

  // sort remaining edges
  std::sort(EI.begin()+l, EI.begin()+l+k, edgeLess);

  // run union-find on remaining edges
  nInMst = unionFindLoop(EI.begin()+l, k, nInMst, UF, mst);

  return mst;
}
