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

#include <math.h>
#include <stddef.h>
#include "parlay/parallel.h"
#include "common/IO.h"
#include "common/geometry.h"
#include "common/geometryIO.h"
#include "geometryData.h"
#include "common/dataGen.h"
#include "common/parseCommandLine.h"
using namespace benchIO;
using namespace dataGen;
using namespace std;

using coord = double;

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"n <outFile>");
  pair<size_t,char*> in = P.sizeAndFileName();
  bool onSphere = P.getOption("-S");

  size_t n = in.first;
  char* fname = in.second;
  parlay::sequence<point3d<coord>> Points(3*n);
  parlay::sequence<tri> Triangles(n);
  coord d = (coord)(1.0/sqrt((double) n));
  parlay::parallel_for (0, n, [&] (size_t i) {
    if (onSphere) Points[3*i] = randOnUnitSphere3d<coord>(i);
    else Points[3*i] = rand3d<coord>(i);
    Points[3*i+1] = Points[3*i] + vector3d<coord>(d,d,0);
    Points[3*i+2] = Points[3*i] + vector3d<coord>(d,0,d);
    Triangles[i][0] = 3*i;
    Triangles[i][1] = 3*i+1;
    Triangles[i][2] = 3*i+2;
  });
  return writeTrianglesToFile(triangles<point3d<coord>>(Points,Triangles),
			      fname);
}
