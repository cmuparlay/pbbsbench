#define NOTMAIN 1
#include <iostream>
#include <limits.h>
#include "graph.h"
#include "parallel.h"
using namespace std;

pbbs::sequence<intT> st(edgeArray<intT> const &EA) {
  edge<intT>* E = EA.E.begin();
  intT m = EA.nonZeros;
  intT n = EA.numRows;
  pbbs::sequence<intT> parents(n, (intT) -1);

  auto find = [&] (intT i) {
    if ((parents[i]) < 0) return i;
    intT j = parents[i];     
    if (parents[j] < 0) return j;
    do j = parents[j]; while (parents[j] >= 0);
    parents[i] = j;
    return j;
  };

  intT* st = pbbs::new_array<intT>(n);
  intT nInSt = 0; 
  for(intT i = 0; i < m; i++){
    intT u = find(E[i].u);
    intT v = find(E[i].v);
    if(u != v){
      if (parents[v] < parents[u]) swap(u,v);
      parents[u] += parents[v];
      parents[v] = u;
      st[nInSt++] = i;
    }
  } 
  return pbbs::sequence<intT>(st, nInSt);
}
