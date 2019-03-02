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
#include "graph.h"
#include "union_find.h"

#include "speculative_for.h"
#include "sequence.h"
#include "parallel.h"
#include "get_time.h"
#include "stlalgs.h"

#include "MST.h"

using namespace std;

// **************************************************************
//    PARALLEL VERSION OF KRUSKAL'S ALGORITHM
// **************************************************************

struct indexedEdge {
  intT u; intT v; intT id;
  indexedEdge(intT u, intT v, intT id) : u(u), v(v), id(id) {}
};

using weight_index = pair<double, intT>;

struct UnionFindStep {
  intT u;  intT v;  
 pbbs::sequence<indexedEdge> const &E;
 pbbs::sequence<reservation> &R;
 unionFind &UF;
 pbbs::sequence<bool> &inST;
  UnionFindStep(pbbs::sequence<indexedEdge> const &E,
		unionFind &UF,
		pbbs::sequence<reservation> &R,
		pbbs::sequence<bool> &inST) :
    E(E), R(R), UF(UF), inST(inST) {}

  bool reserve(intT i) {
    u = UF.find(E[i].u);
    v = UF.find(E[i].v);
    if (u != v) {
      R[v].reserve(i);
      R[u].reserve(i);
      return true;
    } else return false;
  }

  bool commit(intT i) {
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



pbbs::sequence<intT> mst(wghEdgeArray<intT> const &E) { 
  timer t("mst", true);
  size_t m = E.m;
  size_t n = E.n;
  intT k = min<intT>(5 * n / 4, m);
  auto IW = pbbs::delayed_seq<weight_index>(m, [&] (size_t i) {
      return weight_index(E[i].weight, i);});

  auto edgeLess = [&] (weight_index a, weight_index b) { 
    return (a.first == b.first) ? (a.second < b.second) 
    : (a.first < b.first);};
      
  pbbs::random r;
  pbbs::sequence<double> sample(1 + m/1000, [&] (size_t i) -> double {
      return E[r[i]%m].weight;});
  double w = pbbs::sort(sample, std::less<double>())[k/1000];
  t.next("kth smallest");
  
  auto flags = pbbs::delayed_seq<bool>(m, [&] (size_t i) {
      return IW[i].first > w;});
  pbbs::sequence<weight_index> IWS;
  size_t cnt;
  std::tie(IWS, cnt) = pbbs::split_two(IW, flags);
  t.next("partition");

  pbbs::sort_inplace(IWS.slice(0,cnt), edgeLess);
  t.next("first sort");

  pbbs::sequence<reservation> R(n);
  t.next("initialize nodes");

  pbbs::sequence<indexedEdge> Z(cnt, [&] (size_t i) {
      intT j = IWS[i].second;
      return indexedEdge(E[j].u, E[j].v, j);
    });
  t.next("copy to edges");

  pbbs::sequence<bool> mstFlags(m, false);
  unionFind UF(n);
  UnionFindStep UFStep(Z, UF, R,  mstFlags);
  speculative_for(UFStep, 0, cnt, 8);
  t.next("first union find loop");

  auto f = [&] (weight_index wi) {
      intT j = wi.second;
      return UF.find(E[j].u) != UF.find(E[j].v);
  };

  IWS = pbbs::filter(IWS.slice(cnt,m), f);
  t.next("filter out self edges");
  
  pbbs::sort_inplace(IWS.slice(), edgeLess);
  t.next("second sort");

  Z = pbbs::sequence<indexedEdge>(IWS.size(), [&] (size_t i) {
     intT j = IWS[i].second;
     return indexedEdge(E[j].u, E[j].v, j);
    });
  t.next("copy to edges");

  UnionFindStep UFStep2(Z, UF, R, mstFlags);
  speculative_for(UFStep2, 0, IWS.size(), 8);
  t.next("second union find loop");

  pbbs::sequence<intT> mst = pbbs::pack_index<intT>(mstFlags);
  t.next("pack results");

  //cout << "n=" << n << " m=" << m << " nInMst=" << size << endl;
  return mst;
}
