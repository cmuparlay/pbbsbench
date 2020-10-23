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
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/geometry.h"
#include "common/geometryIO.h"
#include "common/parse_command_line.h"
#include "common/get_time.h"
using namespace benchIO;

// *************************************************************
//  SOME DEFINITIONS
// *************************************************************

using coord = double;
using point2 = point2d<coord>;
using point3 = point3d<coord>;

template <class PT, int KK>
struct vertex {
  using pointT = PT;
  int identifier;
  pointT pt;         // the point itself
  vertex* ngh[KK];    // the list of neighbors
  vertex(pointT p, int id) : pt(p), identifier(id) {}
  size_t counter;
};

// *************************************************************
//  TIMING
// *************************************************************

template <int maxK, class point>
void timeNeighbors(parlay::sequence<point> &pts, int k, int rounds, char* outFile) {
  size_t n = pts.size();
  timer t;
  using vtx = vertex<point,maxK>;
  int dimensions = pts[0].dimension();
  auto vv = parlay::tabulate(n, [&] (size_t i) -> vtx {
      return vtx(pts[i],i);
    });
  auto v = parlay::tabulate(n, [&] (size_t i) -> vtx* {
      return &vv[i];});

  // run once for warmup
  ANN<maxK>(v, k);
  for (int i=0; i < rounds; i++) {
    t.start();
    ANN<maxK>(v, k);
    t.next("");
  }

  if (outFile != NULL) {
    int m = n * k;
    parlay::sequence<int> Pout(m);
    parlay::parallel_for (0, n, [&] (size_t i) {
	for (int j=0; j < k; j++)
	  Pout[maxK*i + j] = (v[i]->ngh[j])->identifier;
      });
    writeIntSeqToFile(Pout, outFile);
  }
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-k {1,...,10}] [-d {2,3}] [-o <outFile>] [-r <rounds>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  int rounds = P.getOptionIntValue("-r",1);
  int k = P.getOptionIntValue("-k",1);
  int d = P.getOptionIntValue("-d",2);
  if (k > 10 || k < 1) P.badArgument();
  if (d < 2 || d > 3) P.badArgument();

  if (d == 2) {
    parlay::sequence<point2> PIn = readPointsFromFile<point2>(iFile);
    if (k == 1) timeNeighbors<1>(PIn, 1, rounds, oFile);
    else timeNeighbors<10>(PIn, k, rounds, oFile);
  }

  if (d == 3) {
    parlay::sequence<point3> PIn = readPointsFromFile<point3>(iFile);
    if (k == 1) timeNeighbors<1>(PIn, 1, rounds, oFile);
    else timeNeighbors<10>(PIn, k, rounds, oFile);
  }

}
