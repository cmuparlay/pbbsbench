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
#include <algorithm>
#include <cstring>
#include "parlay/parallel.h"
#include "common/IO.h"
#include "common/graph.h"
#include "common/graphIO.h"
#include "common/parse_command_line.h"
using namespace std;
using namespace benchIO;

// Checks for every vertex if locally maximally matched
int checkMaximalMatching(edgeArray<size_t> const &E, parlay::sequence<size_t> const &EI) {
  size_t m = E.nonZeros;
  size_t nMatches = EI.size();
  size_t n = max(E.numCols,E.numRows);
  parlay::sequence<long> V(n, (long) -1);
  parlay::sequence<bool> flags(m, false);

  parlay::parallel_for (0, nMatches, [&] (size_t i) {
      long idx = EI[i];
      V[E[idx].u] = V[E[idx].v] = idx;
      flags[idx] = 1;
    });

  for (long i=0; i < m; i++) {
    size_t u = E[i].u;
    size_t v = E[i].v;
    if (flags[i]) {
      if (V[u] != i) {
	cout << "maximalMatchingCheck: edges share vertex " << u << endl;
	return 1;
      }
      if (V[v] != i) {
	cout << "maximalMatchingCheck: edges share vertex " << v << endl;
	return 1;
      }
    } else {
      if (u != v && V[u] == -1 && V[v] == -1) {
	cout << "maximalMatchingCheck: neither endpoint matched for edge " << i << endl;
	return 1;
      }
    }
  }
  return 0;
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"<inFile> <outfile>");
  pair<char*,char*> fnames = P.IOFileNames();
  edgeArray<size_t> E = readEdgeArrayFromFile<size_t>(fnames.first);
  parlay::sequence<size_t> Out = readIntSeqFromFile<size_t>(fnames.second);
  return checkMaximalMatching(E, Out);
}
