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
#include "parlay/parallel.h"
#include "common/time_loop.h"
#include "common/graph.h"
#include "common/IO.h"
#include "common/graphIO.h"
#include "common/sequenceIO.h"
#include "common/parse_command_line.h"
#include "BFS.h"
using namespace std;
using namespace benchIO;

void timeBFS(Graph const &G, long source, int rounds, bool verbose, char* outFile) {
  sequence<vertexId> parents;
  time_loop(rounds, 1.0,
	    [&] () {parents.clear();},
	    [&] () {parents = BFS(source, G, verbose);},
	    [&] () {});
  cout << endl;
  if (verbose) {
    size_t visited = parlay::reduce(parlay::delayed_map(parents, [&] (auto p) -> size_t {
       return (p == -1) ? 0 : 1;}));
    cout << "total visited = " << visited << endl;
  }
  if (outFile != NULL) writeSequenceToFile(parents, outFile);
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] [-src source] [-r <rounds>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  int rounds = P.getOptionIntValue("-r",1);
  long source = P.getOptionIntValue("-src",0);
  bool verbose = P.getOption("-v");
  Graph G = readGraphFromFile<vertexId,edgeId>(iFile);
  G.addDegrees();
  timeBFS(G, source, rounds, verbose, oFile);
}
