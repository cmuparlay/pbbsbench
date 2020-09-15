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
#include "parlay/primitives.h"
#include "common/geometry.h"
#include "common/geometryIO.h"
#include "common/topology.h"
#include "common/parse_command_line.h"
#include "refine.h"
#include "common/topology_from_triangles.h"

using std::cout;
using std::endl;
using parlay::sequence;
using parlay::tabulate;
using parlay::reduce;

using namespace benchIO;

#define MIN_ANGLE 30.0

bool skinnyTriangle(triangle<point> *t) {
  if (minAngleCheck(t->vtx[0]->pt, t->vtx[1]->pt, t->vtx[2]->pt, MIN_ANGLE))
    return 1;
  return 0;
}

// double angle(tri *t) {
//   return min(angle(t->vtx[0]->pt, t->vtx[1]->pt, t->vtx[2]->pt),
// 	     min(angle(t->vtx[1]->pt, t->vtx[0]->pt, t->vtx[2]->pt),
// 		 angle(t->vtx[2]->pt, t->vtx[0]->pt, t->vtx[1]->pt)));
// }

bool check(triangles<point> &Tri) {
  size_t m = Tri.numTriangles();
  sequence<vertex<point>> Vertices;
  sequence<triangle<point>> Triangles;
  std::tie(Triangles,Vertices) = topology_from_triangles(Tri);
  if (check_delaunay(Triangles, 10)) return 1;
  
  size_t num_bad = reduce(tabulate(m, [&] (size_t i) -> size_t {
					return skinnyTriangle(&Triangles[i]);}));
  if (num_bad > 0) {
    cout << "Delaunay refine check: " << num_bad << " skinny triangles" << endl;
    return 1;
  }
  return 0;
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,
		"[-r <numtests>] <inFile> <outfile>");
  pair<char*,char*> fnames = P.IOFileNames();
  char* iFile = fnames.first;
  char* oFile = fnames.second;

  triangles<point> T = readTrianglesFromFile<point>(oFile,0);
  return check(T);
}
