#define NOTMAIN 1
#include <iostream>
#include <limits.h>
#include "graph.h"
#include "parallel.h"
using namespace std;

// Assumes root is negative
inline intT find(intT i, intT* parent) {
  if ((parent[i]) < 0) return i;
  intT j = parent[i];     
  if (parent[j] < 0) return j;
  do j = parent[j]; while (parent[j] >= 0);
  parent[i] = j;
  return j;
}

pair<intT*,intT> st(edgeArray<intT> const &EA){
  edge<intT>* E = EA.E.begin();
  intT m = EA.nonZeros;
  intT n = EA.numRows;
  intT *parents = newA(intT,n);
  for(intT i=0;i<n;i++) parents[i] = -1;
  intT* st = newA(intT,m);
  intT nInSt = 0; 
  for(intT i = 0; i < m; i++){
    intT u = find(E[i].u,parents);
    intT v = find(E[i].v,parents);
    if(u != v){
      //union by rank -- join shallower
      //tree to deeper tree
      if(parents[v] < parents[u]) swap(u,v);
      parents[u] += parents[v];
      parents[v] = u;
      st[nInSt++] = i;
    }
  }  
  free(parents); 
  return pair<intT*,intT>(st, nInSt);
}
