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
#include "float.h"
#include <algorithm>
#include <cstring>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/geometry.h"
#include "common/geometryIO.h"
#include "common/parse_command_line.h"
using namespace benchIO;

using coord = double;
using point2 = point2d<coord>;

int checkNeighbors(parlay::sequence<long> &neighbors, parlay::sequence<point2> &P,
		   int k, int r) {
  size_t n = P.size();
  if (neighbors.size() != k * n) {
    std::cout << "error in neighborsCheck: wrong length, n = " << n 
	      << " k = " << k << " neighbors = " << neighbors.size() << std::endl;
    return 1;
  }

  for (long j = 0; j < r; j++) {
    long jj = parlay::hash64(j) % n;

    auto distances = parlay::tabulate(n, [&] (size_t i) -> double {
	return (i == jj) ? DBL_MAX : (P[jj] - P[i]).Length();
      });

    double minD = parlay::reduce(distances, parlay::minm<double>());

    double d = (P[jj] - P[neighbors[k*jj]]).Length();

    double errorTolerance = 1e-6;
    if ((d - minD) / (d + minD)  > errorTolerance) {
      cout << "error in neighborsCheck: for point " << jj 
	   << " min distance reported is: " << d 
	   << " actual is: " << minD << endl;
      return 1;
    }
  }
  return 0;
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,
		"[-k {1,...,100}] [-d {2,3}] [-r <numtests>] <inFile> <outfile>");
  pair<char*,char*> fnames = P.IOFileNames();
  char* iFile = fnames.first;
  char* oFile = fnames.second;

  // number of random points to test
  int r = P.getOptionIntValue("-r",10);

  int k = P.getOptionIntValue("-k",1);
  int d = P.getOptionIntValue("-d",2);
  if (k > 100 || k < 1) P.badArgument();
  if (d < 2 || d > 3) P.badArgument();

  parlay::sequence<long> neighbors = readIntSeqFromFile<long>(oFile);

  if (d == 2) {
    parlay::sequence<point2> PIn = readPointsFromFile<point2>(iFile);
    return checkNeighbors(neighbors, PIn, k, r);
  }
  // } else if (d == 3) {
  //   parlay::sequence<point3d> PIn = readPointsFromFile<point3d>(iFile);
  //   intT n = PIn.n;
  //   point3d* P = PIn.A;
  //   return checkNeighbors(neighbors, PIn.A, PIn.n, k, r);
  // } else return 1;
}
