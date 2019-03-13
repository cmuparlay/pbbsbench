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

// This code will generate (pseudo)-random points in the plummer distribution (3d)
//   or kuzmin distribution (2d)

#include <math.h>
#include "pbbslib/parallel.h"
#include "pbbslib/parse_command_line.h"
#include "common/IO.h"
#include "common/geometry.h"
#include "common/geometryIO.h"
#include "geometryData.h"
using namespace benchIO;
using namespace dataGen;
using namespace std;

using coord = double;

point2d<coord> randKuzmin(size_t i) {
  vector2d<coord> v = vector2d<coord>(randOnUnitSphere2d<coord>(i));
  size_t j = dataGen::hash<size_t>(i);
  double s = dataGen::hash<double>(j);
  double r = sqrt(1.0/((1.0-s)*(1.0-s))-1.0);
  return point2d<coord>(v*r);
}

point3d<coord> randPlummer(size_t i) {
  vector3d<coord> v = vector3d<coord>(randOnUnitSphere3d<coord>(i));
  size_t j = dataGen::hash<size_t>(i);
  double s = pow(dataGen::hash<double>(j),2.0/3.0);
  double r = sqrt(s/(1-s));
  return point3d<coord>(v*r);
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-d {2,3}] n <outFile>");
  pair<size_t,char*> in = P.sizeAndFileName();
  size_t n = in.first;
  char* fname = in.second;
  int dims = P.getOptionIntValue("-d", 2);
  if (dims == 2) {
    pbbs::sequence<point2d<coord>> Points(n, [&] (size_t i) {
	return randKuzmin(i);});
    return writePointsToFile(Points,fname);
  } else if (dims == 3) {
    pbbs::sequence<point3d<coord>> Points(n, [&] (size_t i) {
	return randPlummer(i);});
    return writePointsToFile(Points,fname);
  } 
  P.badArgument();
  return 1;
}
