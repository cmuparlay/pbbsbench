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
#include "ligraLight.h"

using namespace std;

// **************************************************************
//    Using LigraLight
// **************************************************************

parlay::sequence<vertexId> BFS(vertexId start, const Graph &G, bool verbose = false) {
  parlay::internal::timer t("BFS",verbose);
  size_t n = G.numVertices();
  auto parent = parlay::sequence<std::atomic<vertexId>>::from_function(n, [&] (size_t i) {
      return -1;});
  parent[start] = start;

  auto edge_fa = [iparent=parent.begin()] (vertexId u, vertexId v) -> bool {
    vertexId expected = -1;
    return iparent[v].compare_exchange_strong(expected, u);
  };
  auto cond_f = [iparent=parent.begin()] (vertexId v) { return iparent[v] == -1;};
  auto frontier_map = ligra::edge_map(G, edge_fa, cond_f, false, verbose);
  
  auto frontier = ligra::vertex_subset(start);
  t.next("start");

  while (frontier.size() > 0) {
    frontier = frontier_map(frontier);
    t.next("iter");
  }
  return parlay::map(parent, [] (auto const &x) -> vertexId {
      return x.load();});
}
