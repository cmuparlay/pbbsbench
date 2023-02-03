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
#include "../utils/parse_files.h"



#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

using namespace benchIO;

bool report_stats = true;


// *************************************************************
//  TIMING
// *************************************************************
template<typename T>
void timeRange(parlay::sequence<Tvec_point<T>> &pts,
		   parlay::sequence<Tvec_point<T>> &qpoints,
		   int k, int R, int beamSize,
		   int beamSizeQ, double delta, double alpha, double rad, char* outFile,
		   parlay::sequence<ivec_point>& groundTruth, int maxDeg, char* rFile, bool graph_built=false)
{
  size_t n = pts.size();
  auto v = parlay::tabulate(n, [&] (size_t i) -> Tvec_point<T>* {
      return &pts[i];});

  size_t q = qpoints.size();
  auto qpts =  parlay::tabulate(q, [&] (size_t i) -> Tvec_point<T>* {
      return &qpoints[i];});

    time_loop(1, 0,
      [&] () {},
      [&] () {
        RNG<T>(v, k, R, beamSize, beamSizeQ, alpha, delta, rad, qpts, groundTruth, rFile, graph_built);
      },
      [&] () {});

  if(outFile != NULL) {
    std::cout << "Writing graph..."; 
    write_graph(v, outFile, maxDeg); 
    std::cout << " done" << std::endl;
  }

}

// Infile is a file in .fvecs format
int main(int argc, char* argv[]) {
    commandLine P(argc,argv,
    "[-a <alpha>] [-d <delta>] [-R <deg>]"
        "[-L <bm>] [-k <k> ] [-Q <bmq>] [-q <qF>] [-g <gF>] [-o <oF>] [-r <rnds>] [-res [result]] [-b <algoOpt>] [-rad <radius>] <inFile>");

  char* iFile = P.getArgument(0);
  char* qFile = P.getOptionValue("-q");
  char* cFile = P.getOptionValue("-c");
  char* gFile = P.getOptionValue("-g");
  char* oFile = P.getOptionValue("-o");
  char* rFile = P.getOptionValue("-res");
  int R = P.getOptionIntValue("-R", 5);
  if (R < 1) P.badArgument();
  int L = P.getOptionIntValue("-L", 10);
  if (L < 1) P.badArgument();
  int Q = P.getOptionIntValue("-Q", L);
  if (Q < 1) P.badArgument();
  int k = P.getOptionIntValue("-k", 1);
  if (k > 1000 || k < 1) P.badArgument();
  double alpha = P.getOptionDoubleValue("-a", 1.2);
  double delta = P.getOptionDoubleValue("-d", .01);
  int algoOpt = P.getOptionIntValue("-b", 0);
  double radius = P.getOptionDoubleValue("-rad", 96237);

    int maxDeg;
    if(algoOpt == 1) maxDeg = L*R;
    else if(algoOpt == 2) maxDeg = 2*R;
    else maxDeg = R;

    parlay::sequence<ivec_point> groundTruth;
    if (cFile != nullptr)	groundTruth = parse_rangeres(cFile);

    auto [md, points] = parse_uint8bin(iFile, gFile, maxDeg);
    maxDeg = md; 
    auto [fd, qpoints] = parse_uint8bin(qFile, NULL, 0);
    bool graph_built = (gFile != NULL);
    timeRange<uint8_t>(points, qpoints, k, R, L, Q, delta,
          alpha, radius, oFile, groundTruth, maxDeg, rFile, graph_built);

  

}
