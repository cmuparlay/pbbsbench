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
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/geometry.h"
#include "common/geometryIO.h"
#include "common/parse_command_line.h"
#include "common/time_loop.h"
#include "benchUtils.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

using namespace benchIO;


// *************************************************************
//  TIMING
// *************************************************************

template <class point>
void timeNeighbors(parlay::sequence<point> &pts, char* qFile,
  int k, int rounds, int maxDeg, int beamSize, double alpha, char* outFile) 
{
  size_t n = pts.size();
  auto v = parlay::tabulate(n, [&] (size_t i) -> point* {
      return &pts[i];});

  parlay::sequence<fvec_point> qpoints;
  parlay::sequence<fvec_point*> qpts;

  if(qFile != NULL){
    qpoints = parse_fvecs(qFile);
    size_t q = qpoints.size();
    qpts = parlay::tabulate(q, [&] (size_t i) -> point* {
      return &qpoints[i];});

    time_loop(rounds, 1.0,
      [&] () {},
      [&] () {ANN(v, k, maxDeg, beamSize, alpha, qpts);},
      [&] () {});

    if (outFile != NULL) {
      int m = n * (k+1);
      parlay::sequence<int> Pout(m);
      parlay::parallel_for (0, q, [&] (size_t i) {
        Pout[(k+1)*i] = qpts[i]->id;
        for (int j=0; j < k; j++)
          Pout[(k+1)*i + j+1] = (qpts[i]->ngh)[j];
      });
      writeIntSeqToFile(Pout, outFile);
    }
  } else{
      time_loop(rounds, 1.0,
      [&] () {},
      [&] () {ANN(v, k, maxDeg, beamSize, alpha);},
      [&] () {});  

      if (outFile != NULL) {
        int m = n * (k+1);
        parlay::sequence<int> Pout(m);
        parlay::parallel_for (0, n, [&] (size_t i) {
          Pout[(k+1)*i] = v[i]->id;
          for (int j=0; j < k; j++)
            Pout[(k+1)*i + j+1] = (v[i]->ngh)[j];
        });
        writeIntSeqToFile(Pout, outFile);
      }
  }

}

// Infile is a file in .fvecs format
int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-a <alpha>] [-R <maxDeg>] [-L <beamSize>] [-k {1,...,100}] [-o <outFile>] [-q <queryFile>] [-r <rounds>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  char* qFile = P.getOptionValue("-q");
  int rounds = P.getOptionIntValue("-r",1);
  int k = P.getOptionIntValue("-k",1);
  if (k > 100 || k < 1) P.badArgument();
  int maxDeg = P.getOptionIntValue("R", 10);
  int beamSize = P.getOptionIntValue("L", 5);
  double alpha = P.getOptionDoubleValue("alpha", 1.5);

  std::cout << "Input (fvecs format): " << iFile << std::endl;
  auto points = parse_fvecs(iFile);
  // auto qpoints = parse_fvecs(qFile);

  timeNeighbors(points, qFile, k, rounds, maxDeg, beamSize, alpha, oFile);
}
