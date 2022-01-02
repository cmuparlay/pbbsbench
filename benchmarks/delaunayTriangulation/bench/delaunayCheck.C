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
#include "common/topology.h"
#include "common/geometryIO.h"
#include "common/parseCommandLine.h"

#include "delaunay.h"
#include "common/topology_from_triangles.h"
using namespace std;
using namespace benchIO;

using vertex_t = vertex<point>;
using simplex_t = simplex<point>;
using triang_t = triangle<point>;

bool check(triangles<point> &Tri, parlay::sequence<point> &P) {
  size_t m = Tri.numTriangles();
  for (size_t i=0; i < P.size(); i++)
    if (P[i].x != Tri.P[i].x || P[i].y != Tri.P[i].y) {
      cout << "checkDelaunay: prefix of points don't match input at " 
	   << i << endl;
      cout << P[i] << " " << Tri.P[i] << endl;
      return 0;
    }
  auto TriangsVertices = topology_from_triangles(Tri);
  return check_delaunay(TriangsVertices.first, 10);
}
    

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,
		"[-r <numtests>] <inFile> <outfile>");
  pair<char*,char*> fnames = P.IOFileNames();
  char* iFile = fnames.first;
  char* oFile = fnames.second;

  parlay::sequence<point> PIn = readPointsFromFile<point>(iFile);
  triangles<point> T = readTrianglesFromFile<point>(oFile,0);
  check(T, PIn);

  return 0;
}
