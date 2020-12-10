#include <iostream>
#include <limits.h>
#include "parlay/primitives.h"
#include "parlay/parallel.h"
#include "common/get_time.h"
#include "common/graph.h"
#include "common/atomics.h"
#include "algorithm/union_find.h"
#include "ST.h"

parlay::sequence<edgeId> st(edgeArray<vertexId> const &E){
  edgeId m = E.nonZeros;
  vertexId n = E.numRows;
  unionFind<vertexId> UF(n);
  // initialize to an id out of range
  parlay::sequence<edgeId> hooks(n, (edgeId) m);

  parlay::parallel_for (0, m, [&] (edgeId i) {
      vertexId u = E[i].u;
      vertexId v = E[i].v;
      while(1) {
	u = UF.find(u);
	v = UF.find(v);
	if (u == v) break;
	if (u > v) std::swap(u,v);
	if (hooks[u] == m &&
	    pbbs::atomic_compare_and_swap(&hooks[u], m, i)){
	  UF.link(u, v);
	  break;
	}
      }
    }, 1000);

  //get the IDs of the edges in the spanning forest
  parlay::sequence<edgeId> stIdx =  parlay::filter(hooks, [&] (size_t a) {
      return a != m;});
  
  std::cout << "nInSt = " << stIdx.size() << std::endl;
  return stIdx;
}
