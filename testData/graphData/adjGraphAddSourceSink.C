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

// Adds a random integer weight to each node and edge 

#include <math.h>
#include <iostream>

#include "common/graph.h"
#include "common/graphIO.h"
#include "pbbslib/parse_command_line.h"

#include "pbbslib/parallel.h"
using namespace benchIO;
using namespace dataGen;
using namespace std;

typedef long long ll;
typedef wghVertex<intT> V;
ll rand_range(pair<ll,ll> range) {
  return (((ll)rand() << 16) | rand()) % (range.second - range.first + 1) + range.first;
}

int main(int argc, char* argv[]) {
  srand(time(NULL));
  commandLine P(argc,argv,"<inFile> <outFile>");
  pair<char*,char*> fnames = P.IOFileNames();
  char* iFile = fnames.first;
  char* oFile = fnames.second;

  wghGraph<intT> g = readWghGraphFromFile<intT>(iFile);
  intT n = g.n, m = g.m;
  intT maxW = utils::log2Up(n);

  intT S = n, T = n + 1;
  V* v = newA(V, n + 2);
  parallel_for (0, n, [&] (size_t i) {
    v[i] = g.V[i];
    });

  v[S].Neighbors = newA(intT, n);
  v[S].nghWeights = newA(intT, n);
  intT& ds = v[S].degree;
  ds = 0;
  v[T].degree = 0;
  for (intT i = 0; i < n; ++i) {
    intT supply = rand_range(make_pair(-maxW, maxW));
    if (supply > 0) {
      v[S].Neighbors[ds] = i;
      v[S].nghWeights[ds] = supply;
      m++;
      ds++;
    } else if (supply < 0) {
      int d = v[i].degree;
      intT *newneighbors = newA(intT, d + 1);
      intT *neww = newA(intT, d + 1);
      std::copy(v[i].Neighbors, v[i].Neighbors + d, newneighbors);
      std::copy(v[i].nghWeights, v[i].nghWeights + d, neww);
      newneighbors[d] = T;
      neww[d] = -supply;
      v[i].Neighbors = newneighbors;
      v[i].nghWeights = neww;
      v[i].degree++;
      m++;
    }
  }
  wghGraph<intT> gn(v, n + 2, m);
  ofstream out(oFile, ofstream::binary);
  writeFlowGraph<intT>(out, FlowGraph<intT>(gn, S, T));
}
