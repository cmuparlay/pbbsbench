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

#include <limits>
#include "parlay/primitives.h"
#include "parlay/parallel.h"
#include "common/graph.h"
#include "BFS.h"

using namespace std;

// **************************************************************
//    NON DETERMINISTIC BREADTH FIRST SEARCH
// **************************************************************

// **************************************************************
//    Updates the graph so that it is the BFS tree (i.e. the neighbors
//      in the new graph are the children in the bfs tree)
// **************************************************************

std::pair<vertexId,size_t> BFS(vertexId start, Graph &G) {
  vertexId numVertices = G.numVertices();
  edgeId numEdges = G.m;
  vertexId maxIdx = std::numeric_limits<vertexId>::max();

  auto Offsets = parlay::sequence<edgeId>::uninitialized(numVertices+1);
  parlay::sequence<bool> Visited(numVertices, false);
  parlay::sequence<vertexId> Frontier(1, start);

  Visited[start] = true;
  size_t round = 0;
  vertexId totalVisited = 0;
  
  while (Frontier.size() > 0) {
    totalVisited += Frontier.size();
    round++;

    parlay::parallel_for (0, Frontier.size(), [&] (size_t i) {
	Offsets[i] = G[Frontier[i]].degree; });
    
    // Find offsets to write the next frontier for each v in this frontier
    size_t nr = parlay::scan_inplace(Offsets.head(Frontier.size()));
    Offsets[Frontier.size()] = nr;
    auto FrontierNext = parlay::sequence<vertexId>::uninitialized(nr);

    // For each vertex in the frontier try to "hook" unvisited neighbors.
    parlay::parallel_for (0, Frontier.size(), [&] (size_t i) {
	size_t k = 0;
	vertexId v = Frontier[i];
	edgeId o = Offsets[i];
	for (size_t j=0; j < G[v].degree; j++) {
	  vertexId ngh = G[v].Neighbors[j];
	  if (!Visited[ngh] &&  __sync_bool_compare_and_swap(&Visited[ngh], false, true)) {
	    FrontierNext[o+j] = G[v].Neighbors[k++] = ngh;
	  }
	  else FrontierNext[o+j] = -1;
	}
	G.degrees[v] = k;
      });

    // Filter out the empty slots (marked with -1)
    Frontier = parlay::filter(FrontierNext, [&] (vertexId x) {return x >= 0;});
  }

  return std::pair<vertexId,size_t>(totalVisited,round);
}
