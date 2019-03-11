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

#define NOTMAIN 1
#include "sequence.h"
#include "graph.h"
#include "parallel.h"
#include "BFS.h"
#include <limits>
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
  
  pbbs::sequence<edgeId> Offsets(numVertices+1);
  pbbs::sequence<bool> Visited(numVertices, false);
  pbbs::sequence<vertexId> Frontier(1, start);

  Visited[start] = false;
  size_t round = 0;
  vertexId totalVisited = 0;
  
  while (Frontier.size() > 0) {
    totalVisited += Frontier.size();
    round++;

    parallel_for (0, Frontier.size(), [&] (size_t i) {
	Offsets[i] = G[Frontier[i]].degree; });
    
    // Find offsets to write the next frontier for each v in this frontier
    size_t nr = pbbs::scan_inplace(Offsets.slice(0,Frontier.size()), pbbs::addm<edgeId>());
    Offsets[Frontier.size()] = nr;
    pbbs::sequence<vertexId> FrontierNext(nr);

    // For each vertex in the frontier try to "hook" unvisited neighbors.
    parallel_for (0, Frontier.size(), [&] (size_t i) {
	size_t k = 0;
	vertexId v = Frontier[i];
	edgeId o = Offsets[i];
	for (size_t j=0; j < G[v].degree; j++) {
	  vertexId ngh = G[v].Neighbors[j];
	  if (!Visited[ngh] && pbbs::atomic_compare_and_swap(&Visited[ngh], false, true)) {
	    FrontierNext[o+j] = G[v].Neighbors[k++] = ngh;
	  }
	  else FrontierNext[o+j] = -1;
	}
	G.degrees[v] = k;
      });

    // Filter out the empty slots (marked with -1)
    Frontier = pbbs::filter(FrontierNext, [&] (vertexId x) {return x >= 0;});
  }

  return std::pair<vertexId,size_t>(totalVisited,round);
}
