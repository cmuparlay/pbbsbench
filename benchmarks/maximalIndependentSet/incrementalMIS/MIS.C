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
#include "common/speculative_for.h"
#include "common/get_time.h"
#include "MIS.h"
using namespace std;

// **************************************************************
//    MAXIMAL INDEPENDENT SET
// **************************************************************

// For each vertex:
//   Flags = 0 indicates undecided
//   Flags = 1 indicates chosen
//   Flags = 2 indicates a neighbor is chosen
struct MISstep {
  char flag;
  parlay::sequence<char> &Flags;
  Graph const &G;
  MISstep(parlay::sequence<char> & F, const Graph &G) : Flags(F), G(G) {}

  bool reserve(size_t i) {
    size_t d = G[i].degree;
    flag = 1;
    for (size_t j = 0; j < d; j++) {
      size_t ngh = G[i].Neighbors[j];
      if (ngh < i) {
	if (Flags[ngh] == 1) { flag = 2; return 1;}
	// need to wait for higher priority neighbor to decide
	else if (Flags[ngh] == 0) flag = 0; 
      }
    }
    return 1;
  }

  bool commit(size_t i) { return (Flags[i] = flag) > 0;}
};

parlay::sequence<char> maximalIndependentSet(Graph const &GS) {
  size_t n = GS.n;
  parlay::sequence<char> Flags(n, (char) 0);
  MISstep mis(Flags, GS);
  pbbs::speculative_for<vertexId>(mis, 0, n, 20);
  return Flags;
}
