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
#include "common/geometry.h"
#include "common/geometryIO.h"
#include "common/parse_command_line.h"
#include "parlay/primitives.h"
#include "nbody.h"
using namespace std;
using namespace benchIO;

double check(parlay::sequence<point> p, parlay::sequence<vect> forces) {
  size_t n = p.size();
  int nCheck = min<size_t>(n, 200);
  parlay::sequence<double> Err(nCheck);
  double mass = 1.0; // all masses are 1, for now
  
  parlay::parallel_for (0, nCheck, [&] (size_t i) {
    size_t idx = parlay::hash64(i)%n;
    vect force(0.,0.,0.);
    for (size_t j=0; j < n; j++) {
      if (idx != j) {
	vect v = p[j] - p[idx];
	double r2 = v.dot(v);
	force = force + (v * (mass * mass / (r2*sqrt(r2))));
      }
    }
    Err[i] = (force - forces[idx]).Length()/force.Length();
  });
  double total = 0.0;
  for (int i=0; i < nCheck; i++) 
    total += Err[i];
  return total/nCheck;
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv, "[-e <errbound>] <inFile> <outfile>");
  pair<char*,char*> fnames = P.IOFileNames();
  double errorBound = P.getOptionDoubleValue("-e", 1e-6);

  char* iFile = fnames.first;
  char* oFile = fnames.second;

  parlay::sequence<point> PIn = readPointsFromFile<point>(iFile);
  parlay::sequence<vect> POut = readPointsFromFile<vect>(oFile);
  if (PIn.size() != POut.size()) {
    cout << "nbody Check: lengths don't match" << endl;
    return 1;
  }
  double err = check(PIn, POut);
  if (err > errorBound) {
    cout << "nbody Check: RMS error is " << err 
	 << " needs to be less than " << errorBound << endl;
    return 1;
  }
  return 0;
}
