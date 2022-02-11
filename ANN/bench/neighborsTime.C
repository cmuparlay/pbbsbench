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

bool report_stats = true;

// *************************************************************
//  TIMING
// *************************************************************

template<typename T>
void timeNeighbors(parlay::sequence<Tvec_point<T>> &pts,
  int k, int rounds, int maxDeg, int beamSize, double alpha, char* outFile) 
{
  size_t n = pts.size();
  auto v = parlay::tabulate(n, [&] (size_t i) -> Tvec_point<T>* {
      return &pts[i];});

 
  time_loop(rounds, 1.0,
  [&] () {},
  [&] () {ANN<T>(v, k, maxDeg, beamSize, alpha);},
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

template<typename T>
void timeNeighbors(parlay::sequence<Tvec_point<T>> &pts, parlay::sequence<Tvec_point<T>> &qpoints,
  int k, int rounds, int maxDeg, int beamSize, double alpha, char* outFile) 
{
  size_t n = pts.size();
  auto v = parlay::tabulate(n, [&] (size_t i) -> Tvec_point<T>* {
      return &pts[i];});

  size_t q = qpoints.size();
  auto qpts =  parlay::tabulate(q, [&] (size_t i) -> Tvec_point<T>* {   
      return &qpoints[i];});

    time_loop(rounds, 1.0,
      [&] () {},
      [&] () {ANN<T>(v, k, maxDeg, beamSize, alpha, qpts);},
      [&] () {});

    if (outFile != NULL) {
      int m = q * (k+1);
      parlay::sequence<int> Pout(m);
      parlay::parallel_for (0, q, [&] (size_t i) {
        Pout[(k+1)*i] = qpts[i]->id;
        for (int j=0; j < k; j++)
          Pout[(k+1)*i + j+1] = (qpts[i]->ngh)[j];
      });
      writeIntSeqToFile(Pout, outFile);
    }
  

}

// Infile is a file in .fvecs format
int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-a <alpha>] [-R {1,...,1000}] [-L {1,...,1000}] [-k {1,...,100}] [-q <queryFile>] [-o <outFile>] [-r <rounds>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  char* qFile = P.getOptionValue("-q");
  int R = P.getOptionIntValue("-R", 5);
  if (R > 1000 || R < 1) P.badArgument();
  int L = P.getOptionIntValue("-L", 10);
  if (L > 1000 || L < 1) P.badArgument();
  int rounds = P.getOptionIntValue("-r", 1);
  int k = P.getOptionIntValue("-k", 1);
  if (k > 100 || k < 1) P.badArgument();
  double alpha = P.getOptionDoubleValue("-a", 1.5);

  bool fvecs = true;
  std::string filename = std::string(iFile);
  std::string::size_type n = filename.size();
  if(filename[n-5] == 'b') fvecs = false;
  std::cout << "Input (Tvecs format): " << iFile << std::endl;

  if(fvecs){
    parlay::sequence<Tvec_point<float>> points = parse_fvecs(iFile);
    if(qFile != NULL){
      parlay::sequence<Tvec_point<float>> qpoints = parse_fvecs(qFile);
      timeNeighbors<float>(points, qpoints, k, rounds, R, L, alpha, oFile);
    }
    else timeNeighbors<float>(points, k, rounds, R, L, alpha, oFile);
  } 
  else{ 
    parlay::sequence<Tvec_point<uint8_t>> points = parse_bvecs(iFile);
    if(qFile != NULL){
      parlay::sequence<Tvec_point<uint8_t>> qpoints = parse_bvecs(qFile);
      timeNeighbors<uint8_t>(points, qpoints, k, rounds, R, L, alpha, oFile);
    }
    else timeNeighbors<uint8_t>(points, k, rounds, R, L, alpha, oFile);}
}
