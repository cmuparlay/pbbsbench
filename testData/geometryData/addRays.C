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
#include "parlay/primitives.h"
#include "common/IO.h"
#include "common/geometry.h"
#include "common/geometryIO.h"
#include "common/dataGen.h"
#include "common/parse_command_line.h"
#include "geometryData.h"
using namespace benchIO;
using namespace dataGen;
using namespace std;

using point = point3d<float>;
using ray_t = ray<point>;

// p0 is lower corner
// p1 is upper corner
parlay::sequence<ray_t> generateRays(point p0, point p1, size_t n) {
  auto d = p1 - p0;
  auto rays = parlay::tabulate(n, [&] (size_t i) -> ray_t {
    point pl = point(p0.x + d.x * dataGen::hash<double>(4*i+0), 
		     p0.y + d.y * dataGen::hash<double>(4*i+1), 
		     p0.z);
    point pr = point(p0.x + d.x * dataGen::hash<double>(4*i+2), 
		     p0.y + d.y * dataGen::hash<double>(4*i+3), 
		     p1.z);
    return ray_t(pl, pr-pl);
  });
  return rays;
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-n <num>] <triangleInFile> <rayOutFile>");
  pair<char*,char*> fnames = P.IOFileNames();
  char* ifile = fnames.first;
  char* ofile = fnames.second;

  triangles<point> T = readTrianglesFromFile<point>(ifile,1);
  size_t n = P.getOptionLongValue("-n", T.numTriangles());
  auto minf = [] (point a, point b) {return a.minCoords(b);};
  auto maxf = [] (point a, point b) {return a.maxCoords(b);};
  point minPt = parlay::reduce(T.P, parlay::make_monoid(minf,T.P[0]));
  point maxPt = parlay::reduce(T.P, parlay::make_monoid(maxf,T.P[0]));

  // generate as many rays as triangles
  parlay::sequence<ray_t> rays = generateRays(minPt, maxPt, T.numTriangles());
  auto pts = parlay::tabulate(2*n, [&] (size_t i) -> point {
    return (i % 2 == 0) ? rays[i/2].o : point(0,0,0) + rays[i/2].d;});
  return writePointsToFile(pts, ofile);
}
