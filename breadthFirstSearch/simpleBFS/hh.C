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
#include "parlay/internal/get_time.h"
#include "parlay/internal/block_delayed.h"
#include "common/graph.h"
#include "BFS.h"

namespace delayed = parlay::block_delayed;

using namespace std;

// **************************************************************
//    SIMPLE PARALLEL NON DETERMINISTIC BREADTH FIRST SEARCH
// **************************************************************

std::pair<vertexId,size_t> BFS(vertexId start, Graph &G) {
  size_t n = G.numVertices();
  parlay::sequence<std::atomic<vertexId>> parent(n);
  parlay::parallel_for(0, n, [&] (size_t i) {parent[i] = -1;});
  parent[start] = start;

  parlay::sequence<vertexId> frontier(1,start);
  size_t total_visited = 0;
  size_t round = 0;

  while (frontier.size() > 0) {
    total_visited += frontier.size();
    round++;

    auto nested_edges = parlay::map(frontier, [&] (vertexId v) {
	return parlay::delayed_tabulate(G[v].degree, [&, v] (size_t i) {
	    return std::pair(v, G[v].Neighbors[i]);});});
    auto edges = delayed::flatten(nested_edges);

    auto edge_f = [&] (auto u_v) {
      vertexId expected = -1;
      auto [u, v] = u_v;
      return (parent[v] == -1) && parent[v].compare_exchange_strong(expected, u);
    };

    frontier = delayed::filter_map(edges, edge_f, [] (auto x) {return x.second;});
  }
  return std::pair(total_visited, round);
}
