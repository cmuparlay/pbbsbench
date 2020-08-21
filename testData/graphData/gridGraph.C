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
#include "parlay/parallel.h"
#include "common/parse_command_line.h"
#include "common/graph.h"
#include "common/graphIO.h"
#include "common/graphUtils.h"

using namespace benchIO;
using namespace dataGen;
using namespace std;

size_t loc2d(size_t n, size_t i1, size_t i2) {
  return ((i1 + n) % n)*n + (i2 + n) % n;
}

template <class intV>
edgeArray<intV> edge2DMesh(size_t n) {
  size_t dn = round(pow((float) n,1.0/2.0));
  size_t nn = dn*dn;
  size_t nonZeros = 2*nn;
  parlay::sequence<edge<intV>> E(nonZeros);
  parlay::parallel_for (0, dn, [&] (size_t i) {
    for (size_t j=0; j < dn; j++) {
      size_t l = loc2d(dn,i,j);
      E[2*l] = edge<intV>(l,loc2d(dn,i+1,j));
      E[2*l+1] = edge<intV>(l,loc2d(dn,i,j+1));
    }});
  return edgeArray<intV>(std::move(E), nn, nn);
}

size_t loc3d(size_t n, size_t i1, size_t i2, size_t i3) {
  return ((i1 + n) % n)*n*n + ((i2 + n) % n)*n + (i3 + n) % n;
}

template <class intV>
edgeArray<intV> edge3DMesh(size_t n) {
  size_t dn = round(pow((float) n,1.0/3.0));
  size_t nn = dn*dn*dn;
  size_t nonZeros = 3*nn;
  parlay::sequence<edge<intV>> E(nonZeros);
  parlay::parallel_for (0, dn, [&] (size_t i) {
    for (size_t j=0; j < dn; j++) 
      for (size_t k=0; k < dn; k++) {
	size_t l = loc3d(dn,i,j,k);
	E[3*l] =   edge<intV>(l,loc3d(dn,i+1,j,k));
	E[3*l+1] = edge<intV>(l,loc3d(dn,i,j+1,k));
	E[3*l+2] = edge<intV>(l,loc3d(dn,i,j,k+1));
      }});
  return edgeArray<intV>(std::move(E), nn, nn);
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-d {2,3}] [-j] [-o] n <outFile>");
  pair<int,char*> in = P.sizeAndFileName();
  size_t n = in.first;
  char* fname = in.second;
  int dims = P.getOptionIntValue("-d", 2);
  bool ordered = P.getOption("-o");
  bool adjArray = P.getOption("-j");
  edgeArray<size_t> EA;
  if (dims == 2) 
    EA = edge2DMesh<size_t>(n);
  else if (dims == 3) 
    EA = edge3DMesh<size_t>(n);
  else 
    P.badArgument();
  writeGraphFromEdges(EA, fname, adjArray, ordered);
}
