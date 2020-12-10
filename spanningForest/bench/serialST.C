#include <iostream>
#include <limits.h>
#include "common/graph.h"
#include "algorithm/union_find.h"
#include "parlay/primitives.h"
#include "ST.h"
using namespace std;

// vertexId needs to be signed
parlay::sequence<edgeId> st(edgeArray<vertexId> const &E) {
  edgeId m = E.nonZeros;
  vertexId n = E.numRows;
  unionFind<vertexId> UF(n);

  parlay::sequence<edgeId> st(n);
  size_t nInSt = 0; 
  for (edgeId i = 0; i < m; i++){
    vertexId u = UF.find(E[i].u);
    vertexId v = UF.find(E[i].v);
    if (u != v) {
      UF.union_roots(u, v);
      st[nInSt++] = i;
    }
  } 
  return parlay::sequence<edgeId>(st.begin(), st.begin() + nInSt);
}
