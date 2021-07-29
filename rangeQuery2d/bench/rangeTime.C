// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2010 Guy Blelloch and the PBBS team
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
#include "parlay/io.h"
#include "common/time_loop.h"
#include "common/IO.h"
#include "common/geometryIO.h"
#include "common/parse_command_line.h"

#include "range.h"

using namespace std;
using namespace benchIO;

void timeRange(Points const &points, Queries const& queries,
	       int rounds, bool verbose, char* outFile) {
  cout << "start timeRange" << endl;
  long result;
  time_loop(rounds, 2.0,
	    [&] () {},
	    [&] () {result = range(points, queries, verbose);},
	    [&] () {});
  cout << endl;

  cout << "total count = " << result << endl;
  if (outFile != NULL) parlay::chars_to_file(parlay::to_chars(result), outFile);
}

using pointx = point2d<coord>;

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] [-r <rounds>] [-v] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  bool verbose = P.getOption("-v");
  int rounds = P.getOptionIntValue("-r",1);

  parlay::sequence<pointx> A = readPointsFromFile<pointx>(iFile);
  size_t n = A.size();
  size_t num_q = n/3;
  auto points = parlay::map(A.cut(2 * num_q, n), [&] (pointx pt) {return point{pt.x,pt.y};});
  auto queries = parlay::tabulate(num_q, [&] (size_t i) {
     query q;					 
     coord x1 = A[2*i].x;
     coord y1 = A[2*i].y;
     coord x2 = A[2*i+1].x;
     coord y2 = A[2*i+1].y;
     query a{min(x1,x2), max(x1,x2), min(y1,y2), max(y1,y2)};
     return a;});
  timeRange(points, queries, rounds, verbose, oFile);
}

