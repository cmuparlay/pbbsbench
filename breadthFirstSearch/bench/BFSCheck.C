// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011-2019 Guy Blelloch and the PBBS team
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
#include <algorithm>
#include <cstring>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/monoid.h"
#include "common/IO.h"
#include "common/graph.h"
#include "common/graphIO.h"
#include "common/parse_command_line.h"
#include "BFS.h"
using namespace std;
using namespace benchIO;

using vtx_seq = parlay::sequence<vertexId>;

vtx_seq level_number(vtx_seq const &parent) {
  using cp = std::pair<vertexId,vertexId>;
  vertexId n = parent.size();
  auto count_ptr_next = parlay::tabulate(n, [&] (size_t i) -> cp {
      if (parent[i] == -1 || parent[i] == i) return cp(0,i);
      else return cp(1,parent[i]);});
  parlay::sequence<cp> count_ptr;
  // repeated doubling
  while (count_ptr_next != count_ptr) {
    count_ptr = std::move(count_ptr_next);
    count_ptr_next = parlay::map(count_ptr, [&] (cp c_p) {
	auto next = count_ptr[c_p.second];
 	return cp(next.first + c_p.first, next.second);
      });
  }
  return parlay::map(count_ptr_next, [&] (cp c_p) {
      return (parent[c_p.second] == -1) ? n : c_p.first;});
}

// Checks if parents is valid BFS tree relative to G starting at start
int checkBFS(size_t start, vtx_seq const &parents, Graph const &G) {
  size_t n = G.numVertices();
  if (n != parents.size()) {
    cout << "BFSCheck: vertex counts don't match: " << n << ", " << parents.size() << endl;
    return 1;
  }

  parlay::sequence level = level_number(parents);
  auto min_ngh_level = [&] (vertexId v) {
    if (G[v].degree == 0) return (vertexId) level.size();
    auto nghs = parlay::make_slice(G[v].Neighbors, G[v].Neighbors+G[v].degree);
    auto ngh_level = parlay::delayed_map(nghs, [&] (vertexId u) {
	return level[u];});
    return parlay::reduce(ngh_level, parlay::minm<vertexId>());
  };

  auto flags = parlay::tabulate(n, [&] (vertexId v) {
      if (v == start) return (parents[v] == v);
      if (v == start) return (parents[v] == v);
      auto ngh_level = min_ngh_level(v);
      return (level[v] == n && ngh_level == n ||
	      level[v] == ngh_level + 1);
    });
  auto err = parlay::find(flags, false) - flags.begin();
  if (err < n) {
    cout << "BFSCheck: failed at vertex " << err
	 << ": level is " << level[err]
	 << ", min neighbor level is " << min_ngh_level(err) << endl;
    return 1;
  }
  return 0;
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"<inFile> <outfile>");
  pair<char*,char*> fnames = P.IOFileNames();
  char* iFile = fnames.first;
  char* oFile = fnames.second;

  Graph G = readGraphFromFile<vertexId,edgeId>(iFile);
  vtx_seq parents = readIntSeqFromFile<vertexId>(oFile);

  return checkBFS(0, parents, G);
}
