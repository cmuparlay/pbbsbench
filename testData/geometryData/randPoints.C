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

// This code will generate uniformly (pseudo)-random points
//   By default they will be in the unit cube [-1..1] in each dimension
//   The -s argument will place them in a unit sphere centered at 0 with
//      unit radius
//   The -S argument will place them on the surface of the unit sphere
//   Only one of -s or -S should be used

#include <math.h>
#include "parlay/parallel.h"
#include "common/parse_command_line.h"
#include "common/IO.h"
#include "common/geometry.h"
#include "common/geometryIO.h"
#include "geometryData.h"
using namespace benchIO;
using namespace dataGen;
using namespace std;

using coord = double;

template <class coord>
point2d<coord> randKuzmin(size_t i) {
  vector2d<coord> v = vector2d<coord>(randOnUnitSphere2d<coord>(i));
  size_t j = dataGen::hash<size_t>(i);
  double s = dataGen::hash<double>(j);
  double r = sqrt(1.0/((1.0-s)*(1.0-s))-1.0);
  return point2d<coord>(v*r);
}

template <class coord>
point3d<coord> randPlummer(size_t i) {
  vector3d<coord> v = vector3d<coord>(randOnUnitSphere3d<coord>(i));
  size_t j = dataGen::hash<size_t>(i);
  double s = pow(dataGen::hash<double>(j),2.0/3.0);
  double r = sqrt(s/(1-s));
  return point3d<coord>(v*r);
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-s] [-S] [-k] [-p] [-d {2,3}] n <outFile>\n");
  pair<size_t, char*> in = P.sizeAndFileName();
  size_t n = in.first;
  char* fname = in.second;
  int dims = P.getOptionIntValue("-d", 2);
  bool inSphere = P.getOption("-s");
  bool onSphere = P.getOption("-S");
  bool plummerOrKuzmin = P.getOption("-k") || P.getOption("-p");

  if (dims == 2) {
    auto Points = parlay::tabulate(n, [&] (size_t i) -> point2d<coord> {
	if (inSphere) return randInUnitSphere2d<coord>(i);
	else if (onSphere) return randOnUnitSphere2d<coord>(i);
	else if (plummerOrKuzmin) return randKuzmin<coord>(i);
	else return rand2d<coord>(i);
      });
    return writePointsToFile(Points,fname);
  } else if (dims == 3) {
    auto Points = parlay::tabulate(n, [&] (size_t i) -> point3d<coord> {
	if (inSphere) return randInUnitSphere3d<coord>(i);
	else if (onSphere) return randOnUnitSphere3d<coord>(i);
	else if (plummerOrKuzmin) return randPlummer<coord>(i);
	else return rand3d<coord>(i);
      });
    return writePointsToFile(Points,fname);
  } 
  P.badArgument();
  return 1;
}
