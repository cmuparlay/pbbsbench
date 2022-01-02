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
#include "common/geometry.h"
#include "common/geometryIO.h"
#include "common/parse_command_line.h"
using namespace std;
using namespace benchIO;

using coord = double;
using point = point2d<coord>;

bool checkHull(parlay::sequence<point> const &P, parlay::sequence<size_t> const &I) {
  size_t n = P.size();
  size_t nOut = I.size();
  auto PO = parlay::tabulate(nOut, [&] (size_t i) -> point {return P[I[i]];});
  auto pless = [&] (point a, point b) {
    return (a.x < b.x) || ((a.x == b.x) && (a.y < b.y));};
  auto eq = [&] (point a, point b) {return (a.x == b.x) && (a.y == b.y);};
  size_t idx = parlay::max_element(PO, pless) - PO.begin();
  parlay::sequence<point> PS= parlay::sort(P, pless);
  if (!eq(PS[0], PO[0])) { 
    cout << "checkHull: bad leftmost point" << endl;
    PS[0].print();  PO[0].print(); cout << endl;
    return 1;
  }
  if (!eq(PS[n-1], PO[idx])) {
    cout << "checkHull: bad rightmost point" << endl;
    return 1;
  }
  size_t k = 1;
  for (size_t i = 0; i < idx; i++) {
    if (i > 0 && counterClockwise(PO[i-1],PO[i],PS[i+1])) {
	cout << "checkHull: not convex" << endl;
	return 1;
    }
    if (PO[i].x > PO[i+1].x) {
      cout << "checkHull: not sorted by x" << endl;
      return 1;
    }
    while (!eq(PS[k], PO[i+1]) && k < n)
      if (counterClockwise(PO[i],PO[i+1],PS[k++])) {
	cout << "checkHull: above hull" << endl;
	return 1;
      }
    if (k == n) {
      cout << "checkHull: unexpected points in hull" << endl;
      return 1;
    }
    k++;
  }
  return 0;
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"<inFile> <outfile>");
  pair<char*,char*> fnames = P.IOFileNames();
  char* iFile = fnames.first;
  char* oFile = fnames.second;

  parlay::sequence<point> PIn = readPointsFromFile<point>(iFile);
  parlay::sequence<size_t> I = readIntSeqFromFile<size_t>(oFile);
  return checkHull(PIn, I);
}
