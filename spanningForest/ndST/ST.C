#define NOTMAIN 1
#include <iostream>
#include <limits.h>
#include "sequence.h"
#include "get_time.h"
#include "graph.h"
#include "parallel.h"
#include "union_find.h"
#include "ST.h"
using namespace std;

constexpr vertexId max_vertexId = std::numeric_limits<vertexId>::max();

pbbs::sequence<vertexId> st(edgeArray<vertexId> const &E){
  size_t m = E.nonZeros;
  size_t n = E.numRows;
  pbbs::sequence<vertexId> parents(n, (vertexId) -1);
  pbbs::sequence<vertexId> hooks(n, max_vertexId);

  // Assumes root is negative 
  // Not making parent array volatile improves
  // performance and doesn't affect correctness
  auto find = [&] (vertexId i) {
    vertexId j = i; 
    if (parents[j] < 0) return j;
    do j = parents[j];
    while (parents[j] >= 0);
    //note: path compression can happen in parallel in the same tree, so
    //only link from smaller to larger to avoid cycles
    vertexId tmp;
    while((tmp = parents[i]) <j ){
      parents[i] = j; i = tmp;} 
    return j;
  };
  
  parallel_for (0, m, [&] (vertexId i) {
      vertexId u = E[i].u, v = E[i].v;
      while(1) {
	u = find(u);
	v = find(v);
	if(u == v) break;
	if(u > v) swap(u,v);
	//if successful, store the ID of the edge used in hooks[u]
	if(hooks[u] == max_vertexId &&
	   pbbs::atomic_compare_and_swap(&hooks[u], max_vertexId, i)){
	  parents[u]=v;
	  break;
	}
      }
    }, 1000);

  //get the IDs of the edges in the spanning forest
  pbbs::sequence<vertexId> stIdx =  pbbs::filter(hooks, [&] (size_t a) {
      return a != max_vertexId;});
  
  cout<< "nInSt = " << stIdx.size() << endl;
  return stIdx;
}
