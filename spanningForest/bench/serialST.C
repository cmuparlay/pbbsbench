#define NOTMAIN 1
#include <iostream>
#include <limits.h>
#include "graph.h"
#include "parallel.h"
#include "ST.h"
#include "union_find.h"
using namespace std;

// vertexId needs to be signed
pbbs::sequence<edgeId> st(edgeArray<vertexId> const &E) {
  edgeId m = E.nonZeros;
  vertexId n = E.numRows;
  unionFind<vertexId> UF(n);

  edgeId* st = pbbs::new_array<edgeId>(n);
  size_t nInSt = 0; 
  for (edgeId i = 0; i < m; i++){
    vertexId u = UF.find(E[i].u);
    vertexId v = UF.find(E[i].v);
    if (u != v) {
      UF.union_roots(u, v);
      st[nInSt++] = i;
    }
  } 
  return pbbs::sequence<edgeId>(st, nInSt);
}
