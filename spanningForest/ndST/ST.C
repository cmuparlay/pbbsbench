#define NOTMAIN 1
#include <iostream>
#include <limits.h>
#include "sequence.h"
#include "get_time.h"
#include "graph.h"
#include "parallel.h"
#include "union_find.h"
using namespace std;

constexpr intT max_intT = std::numeric_limits<intT>::max();

pair<intT*, intT> st(edgeArray<intT> const &E){
  intT m = E.nonZeros;
  intT n = E.numRows;
  pbbs::sequence<intT> parents(n, (intT) -1);
  pbbs::sequence<intT> hooks(n, max_intT);

  // Assumes root is negative 
  // Not making parent array volatile improves
  // performance and doesn't affect correctness
  auto find = [&] (intT i) {
    intT j = i; 
    if (parents[j] < 0) return j;
    do j = parents[j];
    while (parents[j] >= 0);
    //note: path compression can happen in parallel in the same tree, so
    //only link from smaller to larger to avoid cycles
    intT tmp;
    while((tmp = parents[i]) <j ){
      parents[i] = j; i = tmp;} 
    return j;
  };
  
  parallel_for (0, m, [&] (intT i) {
      intT u = E[i].u, v = E[i].v;
      while(1) {
	u = find(u);
	v = find(v);
	if(u == v) break;
	if(u > v) swap(u,v);
	//if successful, store the ID of the edge used in hooks[u]
	if(hooks[u] == max_intT &&
	   pbbs::atomic_compare_and_swap(&hooks[u], max_intT, i)){
	  parents[u]=v;
	  break;
	}
      }
    }, 1000);

  //get the IDs of the edges in the spanning forest
  pbbs::sequence<intT> stIdx =  pbbs::filter(hooks, [&] (intT a) {
      return a != max_intT;});
  
  size_t outSize = stIdx.size();
  cout<< "nInSt = " << outSize << endl;
  return pair<intT*,intT>(stIdx.to_array(), outSize);
}
