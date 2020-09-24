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
#include "parlay/primitives.h"
#include "common/geometry.h"
#include "common/geometryIO.h"
#include "common/parse_command_line.h"
#include "ray.h"
using namespace std;
using namespace benchIO;

// WARNING : THIS IS CURRENTLY JUST A STUB

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"<triangleFile> <rayFile> <intersectFile>");
  char* triFile = P.getArgument(2);
  char* rayFile = P.getArgument(1);
  char* intersectFile = P.getArgument(0);

  // the 1 argument means that the vertices are labeled starting at 1
  triangles<point> T = readTrianglesFromFile<point>(triFile,1);
  parlay::sequence<point> Pts = readPointsFromFile<point>(rayFile);
  parlay::sequence<size_t> Is = readIntSeqFromFile<size_t>(intersectFile);

  return 0;
}
