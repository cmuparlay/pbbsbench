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
#include "common/time_loop.h"
#include "common/geometry.h"
#include "common/geometryIO.h"
#include "common/parse_command_line.h"
#include "parlay/primitives.h"
#include "ray.h"

using namespace std;
using namespace benchIO;

void timeRayCast(triangles<point> T, parlay::sequence<ray<point>> rays, 
		 int rounds, char* outFile) {
  parlay::sequence<index_t> R;
  time_loop(rounds, 2.0,
	    [&] () {R.clear();},
	    [&] () {R = rayCast(T, rays, false);},
	    [&] () {});
  cout << endl;
  if (outFile != NULL) writeIntSeqToFile(R, outFile);
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] [-r <rounds>] <triangleFile> <rayFile>");
   pair<char*,char*> fnames = P.IOFileNames();
  char* triFile = fnames.first;
  char* rayFile = fnames.second;
  char* oFile = P.getOptionValue("-o");
  int rounds = P.getOptionIntValue("-r",1);

  // the 1 argument means that the vertices are labeled starting at 1
  triangles<point> T = readTrianglesFromFile<point>(triFile, 1);
  parlay::sequence<point> Pts = readPointsFromFile<point>(rayFile);
  size_t n = Pts.size()/2;
  auto rays = parlay::tabulate(n, [&] (size_t i) -> ray<point> {
      return ray<point>(Pts[2*i], Pts[2*i+1]-point(0,0,0));});
  timeRayCast(T, rays, rounds, oFile);
}
