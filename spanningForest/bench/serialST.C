#define NOTMAIN 1
#include <iostream>
#include <limits.h>
#include "graph.h"
#include "parallel.h"
#include "ST.h"
using namespace std;

pbbs::sequence<vertexId> st(edgeArray<vertexId> const &EA) {
  edge<vertexId>* E = EA.E.begin();
  size_t m = EA.nonZeros;
  size_t n = EA.numRows;
  pbbs::sequence<vertexId> parents(n, (vertexId) -1);

  auto find = [&] (vertexId i) {
    if ((parents[i]) < 0) return i;
    vertexId j = parents[i];     
    if (parents[j] < 0) return j;
    do j = parents[j]; while (parents[j] >= 0);
    parents[i] = j;
    return j;
  };

  vertexId* st = pbbs::new_array<vertexId>(n);
  size_t nInSt = 0; 
  for (size_t i = 0; i < m; i++){
    vertexId u = find(E[i].u);
    vertexId v = find(E[i].v);
    if(u != v){
      if (parents[v] < parents[u]) swap(u,v);
      parents[u] += parents[v];
      parents[v] = u;
      st[nInSt++] = i;
    }
  } 
  return pbbs::sequence<vertexId>(st, nInSt);
}
