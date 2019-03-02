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
#include <algorithm>
#include "graph.h"
#include "get_time.h"
#include "MST.h"
#include "parallel.h"
using namespace std;


// **************************************************************
//    FIND OPERATION FOR UNION FIND
// **************************************************************

// Assumes root is negative
intT find(intT i, intT* parent) {
  if (parent[i] < 0) return i;
  else {
    intT j = i;
    intT tmp;

    // Find root
    do j = parent[j]; while (parent[j] >= 0);

    // path compress
    while ((tmp = parent[i]) != j) { parent[i] = j; i = tmp;}

    return j;
  }
}

// **************************************************************
//    SERIAL MST USING KRUSKAL's ALGORITHM WITH FILTERING
// **************************************************************

struct edgei {
  intT u;
  intT v;
  double weight;
  intT id;
  edgei() {}
  edgei(intT _u, intT _v, double w, intT _id) 
    : u(_u), v(_v), id(_id), weight(w) {}
};

struct edgeLess : std::binary_function <edgei ,edgei ,bool> {
  bool operator() (edgei const& a, edgei const& b) const
  { return (a.weight == b.weight) ? (a.id < b.id) : (a.weight < b.weight); }
};

int unionFindLoop(edgei* E, intT m, intT nInMst, intT* parent, intT* mst) {
  for (intT i = 0; i < m; i++) {
    intT u = find(E[i].u, parent);
    intT v = find(E[i].v, parent);
    
    // union operation 
    if (u != v) {
      if (parent[v] < parent[u]) swap(u,v);
      parent[u] += parent[v]; 
      parent[v] = u;
      mst[nInMst++] = E[i].id;
    }
  }
  return nInMst;
}

pbbs::sequence<intT> mst(wghEdgeArray<intT> const &G) { 
  edgei* EI = pbbs::new_array<edgei>(G.m);
  for (intT i=0; i < G.m; i++) 
    EI[i] = edgei(G.E[i].u, G.E[i].v, G.E[i].weight, i);

  intT l = min(4*G.n/3,G.m);
  std::nth_element(EI, EI+l, EI+G.m, edgeLess());

  std::sort(EI, EI+l, edgeLess());

  intT *parent = pbbs::new_array<intT>(G.n);
  for (intT i=0; i < G.n; i++) parent[i] = -1;

  intT *mst = pbbs::new_array<intT>(G.n);
  intT nInMst = unionFindLoop(EI, l, 0, parent, mst);

  edgei *f = EI+l;
  for (edgei*e = EI+l; e < EI + G.m; e++) {
    intT u = find(e->u, parent);
    intT v = find(e->v, parent);
    if (u != v) *f++ = *e;
  }
  intT k = f - (EI+l);

  std::sort(EI+l, f, edgeLess());

  nInMst = unionFindLoop(EI+l, k, nInMst, parent, mst);

  //cout << "n=" << G.n << " m=" << G.m << " nInMst=" << nInMst << endl;
  pbbs::free_array(EI);
  pbbs::free_array(parent);
  return pbbs::sequence<intT>(mst, nInMst);
}
