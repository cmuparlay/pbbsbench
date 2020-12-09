// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011-2019 Guy Blelloch and the PBBS team
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
#include "parlay/internal/get_time.h"
#include "common/graph.h"
#include "common/speculative_for.h"
#include "algorithm/kth_smallest.h"
#include "algorithm/union_find.h"
#include "MST.h"

using namespace std;

// **************************************************************
//    PARALLEL VERSION OF KRUSKAL'S ALGORITHM
// **************************************************************

// need to tag each edge with an index so we can keep track of
// which edges are added to the MST
struct indexedEdge {
  vertexId u; vertexId v; edgeId id; edgeWeight w;
  indexedEdge(vertexId u, vertexId v, edgeId id, edgeWeight w)
    : u(u), v(v), id(id), w(w){}
  indexedEdge() {};
};

using reservation = pbbs::reservation<edgeId>;

struct UnionFindStep {
  parlay::sequence<indexedEdge> &E;
  parlay::sequence<reservation> &R;
  unionFind<vertexId> &UF;
  parlay::sequence<bool> &inST;
  UnionFindStep(parlay::sequence<indexedEdge> &E,
		unionFind<vertexId> &UF,
		parlay::sequence<reservation> &R,
		parlay::sequence<bool> &inST) :
    E(E), R(R), UF(UF), inST(inST) {}

  bool reserve(edgeId i) {
    vertexId u = E[i].u = UF.find(E[i].u);
    vertexId v = E[i].v = UF.find(E[i].v);
    if (u != v) {
      R[v].reserve(i);
      R[u].reserve(i);
      return true;
    } else return false;
  }

  bool commit(edgeId i) {
    vertexId u = E[i].u;
    vertexId v = E[i].v;
    if (R[v].check(i)) {
      R[u].checkReset(i); 
      UF.link(v, u); 
      inST[E[i].id] = true;
      return true;}
    else if (R[u].check(i)) {
      UF.link(u, v);
      inST[E[i].id] = true;
      return true; }
    else return false;
  }
};

parlay::sequence<edgeId> mst(wghEdgeArray<vertexId,edgeWeight> &E) { 
  parlay::internal::timer t("mst", true);
  size_t m = E.m;
  size_t n = E.n;
  size_t k = min<size_t>(5 * n / 4, m);

  // equal edge weights will prioritize the earliest one
  auto edgeLess = [&] (indexedEdge a, indexedEdge b) { 
    return (a.w < b.w) || ((a.w == b.w) && (a.id < b.id));};

  // tag each edge with an index
  auto IW = parlay::delayed_seq<indexedEdge>(m, [&] (size_t i) {
      return indexedEdge(E[i].u, E[i].v, i, E[i].weight);});

  indexedEdge kth = pbbs::approximate_kth_smallest(IW, k, edgeLess);
  t.next("approximate kth smallest");
  
  auto IW1 = parlay::filter(IW, [&] (indexedEdge e) {
      return edgeLess(e, kth);}); //edgeLess(e, kth);});
  t.next("filter those less than kth smallest as prefix");
  
  parlay::sort_inplace(IW1, edgeLess);
  t.next("sort prefix");

  parlay::sequence<bool> mstFlags(m, false);
  unionFind<vertexId> UF(n);
  parlay::sequence<reservation> R(n);
  UnionFindStep UFStep1(IW1, UF, R,  mstFlags);
  pbbs::speculative_for<vertexId>(UFStep1, 0, IW1.size(), 5, false);
  t.next("union find loop on prefix");

  auto IW2 = parlay::filter(IW, [&] (indexedEdge e) {
      return !edgeLess(e, kth) && UF.find(e.u) != UF.find(e.v);});
  t.next("filter those that are not self edges");
  
  parlay::sort_inplace(IW2, edgeLess);
  t.next("sort remaining");

  UnionFindStep UFStep2(IW2, UF, R, mstFlags);
  pbbs::speculative_for<vertexId>(UFStep2, 0, IW2.size(), 5, false);
  t.next("union find loop on remaining");

  parlay::sequence<edgeId> mst = parlay::internal::pack_index<edgeId>(mstFlags);
  t.next("pack out results");

  //cout << "n=" << n << " m=" << m << " nInMst=" << size << endl;
  return mst;
}

