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

#include <iostream>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/graph.h"
#include "common/atomics.h"
#include "MIS.h"
using namespace std;

// **************************************************************
//    MAXIMAL INDEPENDENT SET
// **************************************************************

// Flags[v] = 0 means undecided, Flags[v] = 1 means vertex is in MIS,
// and Flags[v] = 2 means vertex is not in MIS

// nondeterministic greed algorithm
parlay::sequence<char> maximalIndependentSet(const Graph &G) {
  size_t n = G.n;
  parlay::sequence<char> Flags(n, (char) 0);
  parlay::sequence<bool> V(n, false);
  
  parlay::parallel_for(0, n, [&] (size_t i) {
      size_t v = i;
      while (1) {
	//drop out if already in or out of MIS
	if (Flags[v]) break;
	//try to lock self and neighbors
	if (pbbs::atomic_compare_and_swap<bool>(&V[v], false, true)) {
	  size_t k = 0;
	  for (size_t j = 0; j < G[v].degree; j++){
	    vertexId ngh = G[v].Neighbors[j];
	    // if ngh is not in MIS or we successfully 
	    // acquire lock, increment k
	    if (Flags[ngh] == 2 || pbbs::atomic_compare_and_swap(&V[ngh], false, true))
	      k++;
	    else break;
	  }
	  if(k == G[v].degree){ 
	    //win on self and neighbors so fill flags
	    Flags[v] = 1;
	    for(size_t j = 0; j < G[v].degree; j++){
	      vertexId ngh = G[v].Neighbors[j]; 
	      if(Flags[ngh] != 2) Flags[ngh] = 2;
	    }
	  } else { 
	    //lose so reset V values up to point
	    //where it lost
	    V[v] = false;
	    for(size_t j = 0; j < k; j++){
	      vertexId ngh = G[v].Neighbors[j];
	      if(Flags[ngh] != 2) V[ngh] = false;
	    }
	  }
	}
      }
    });
  return Flags;
}
