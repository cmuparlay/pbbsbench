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
#include "parlay/primitives.h"
#include "common/IO.h"
#include "common/graph.h"
#include "common/graphIO.h"
#include "common/parse_command_line.h"

using namespace std;
using namespace benchIO;

// Checks if valid maximal independent set
int checkMaximalIndependentSet(graph<size_t> const &G,
			       parlay::sequence<size_t> const &Flags) {
  size_t n = G.n;
  using P = std::pair<size_t,size_t>;
  P R = parlay::reduce(parlay::delayed_seq<P>(n, [&] (size_t i) -> P {
	bool hasNeighbor = false;
	for (size_t j=0; j < G[i].degree; j++) {
	  size_t ngh = G[i].Neighbors[j];
	  if (Flags[ngh] == 1) {
	    hasNeighbor = true;
	    if (Flags[i] == 1) return P(i, ngh);
	  }
	}
	if ((Flags[i] != 1) && !hasNeighbor) return P(i,n);
	return P(n,n);
      }), parlay::minm<P>());
  if (R.first < n) {
    if (R.second < n) 
      cout << "checkMaximalIndependentSet: bad edge ("
	   << R.first << ", " << R.second << ")" << endl;
    else
      cout << "checkMaximalIndependentSet: bad vertex "
	   << R.first << endl;
    return 1;
  }      
  return 0;
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"<inFile> <outfile>");
  pair<char*,char*> fnames = P.IOFileNames();
  char* iFile = fnames.first;
  char* oFile = fnames.second;

  graph<size_t> G = readGraphFromFile<size_t>(iFile);
  parlay::sequence<size_t> Out = readIntSeqFromFile<size_t>(oFile);
  if (Out.size() != G.n) {
    cout << "checkMaximalIndependentSet: output file not of right length" << endl;
    return(1);
  }

  return checkMaximalIndependentSet(G, Out);
}
