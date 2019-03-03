#define NOTMAIN 1
#include <iostream>
#include <limits.h>
#include "sequence.h"
#include "get_time.h"
#include "graph.h"
#include "parallel.h"
#include "union_find.h"
#include "ST.h"

pbbs::sequence<edgeId> st(edgeArray<vertexId> const &E){
  edgeId m = E.nonZeros;
  vertexId n = E.numRows;
  unionFind<vertexId> UF(n);
  // initialize to an id out of range
  pbbs::sequence<edgeId> hooks(n, (edgeId) m);

  parallel_for (0, m, [&] (edgeId i) {
      vertexId u = E[i].u;
      vertexId v = E[i].v;
      while(1) {
	u = UF.find(u);
	v = UF.find(v);
	if (u == v) break;
	if (u > v) std::swap(u,v);
	//if successful, store the ID of the edge used in hooks[u]
	if (hooks[u] == m &&
	    pbbs::atomic_compare_and_swap(&hooks[u], m, i)){
	  UF.link(u, v);
	  break;
	}
      }
    }, 1000);

  //get the IDs of the edges in the spanning forest
  pbbs::sequence<edgeId> stIdx =  pbbs::filter(hooks, [&] (size_t a) {
      return a != m;});
  
  cout<< "nInSt = " << stIdx.size() << endl;
  return stIdx;
}
