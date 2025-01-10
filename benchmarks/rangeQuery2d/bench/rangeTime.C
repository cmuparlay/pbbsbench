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
#include "common/sequenceIO.h"
#include "common/geometryIO.h"
#include "common/parse_command_line.h"

#include "range.h"

using namespace std;
using namespace benchIO;
using parlay::sequence;
void timeRange(Points const &points, Queries const& queries,
	       int rounds, bool verbose, char* outFile) {
  cout << "start timeRange" << endl;

  if (queries.size()< 100){
    parlay::map(queries,[&] (query q) {printf("query: (%f, %f) (%f, %f)\n",q.x1,q.y1,q.x2,q.y2 ); return 0;});
  }
  if (points.size()< 100){
    parlay::map(points,[&] (point p) {printf("point: (%f, %f) \n",p.x,p.y ); return 0;});
  }
  sequence<long> R;
  time_loop(rounds, 2.0,
	    [&] () {R.clear();},
	    [&] () {R = range(points, queries, verbose);},
	    [&] () {});
  cout << endl;
  long total = parlay::reduce(R);
  cout << "total count = " << total << endl;
  if (outFile != NULL) writeSequenceToFile(R, outFile);
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

