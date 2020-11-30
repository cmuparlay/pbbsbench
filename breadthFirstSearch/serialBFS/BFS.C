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

#include "common/graph.h"
#include "parlay/sequence.h"
#include "BFS.h"
using namespace std;

// **************************************************************
//    SERIAL BREADTH FIRST SEARCH
// **************************************************************

parlay::sequence<vertexId> BFS(vertexId start, const Graph &G,
			       bool verbose = false) {
  size_t n = G.numVertices();
  parlay::sequence<vertexId> Frontier(n);
  parlay::sequence<vertexId> Parent(n, -1);

  size_t bot = 0;
  size_t top = 1;
  Frontier[0] = start;
  Parent[start] = start;

  while (top > bot) {
    vertexId v = Frontier[bot++];
    size_t k = 0;
    for (size_t j = 0; j < G[v].degree; j++) {
      vertexId ngh = G[v].Neighbors[j];
      if (Parent[ngh] == -1) {
	Frontier[top++] = ngh;
	Parent[ngh] = v;
      }
    }
  }
  return Parent;
}


