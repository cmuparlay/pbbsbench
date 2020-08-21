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

#include "parlay/parallel.h"
#include "common/parse_command_line.h"
#include "common/get_time.h"
#include "common/graph.h"
#include "common/graphIO.h"
#include "common/graphUtils.h"
using namespace benchIO;
using namespace dataGen;
using namespace std;

// Generates an undirected graph with n vertices with approximately degree 
// neighbors per vertex.
// Edges  are distributed so they appear to come from
// a dim-dimensional space.   In particular an edge (i,j) will have
// probability roughly proportional to (1/|i-j|)^{(d+1)/d}, giving 
// separators of size about n^{(d-1)/d}.    
template <class intV>
edgeArray<intV> edgeRandomWithDimension(size_t dim, size_t degree, size_t numRows) {
  size_t nonZeros = numRows*degree;
  auto E = parlay::tabulate(nonZeros, [&] (size_t k) -> edge<intV> {
      size_t i = k / degree;
      size_t j;
      if (dim==0) {
	size_t h = k;
	do {
	  j = ((h = dataGen::hash<intV>(h)) % numRows);
	} while (j == i);
      } else {
	size_t pow = dim+2;
	size_t h = k;
	do {
	  while ((((h = dataGen::hash<intV>(h)) % 1000003) < 500001)) pow += dim;
	  j = (i + ((h = dataGen::hash<intV>(h)) % (((long) 1) << pow))) % numRows;
	} while (j == i);
      }
      return edge<intV>(i, j);
    });
  return edgeArray<intV>(std::move(E), numRows, numRows);
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-m <numedges>] [-d <dims>] [-o] [-j] n <outFile>");
  pair<size_t,char*> in = P.sizeAndFileName();
  size_t n = in.first;
  char* fname = in.second;
  int dim = P.getOptionIntValue("-d", 0);
  size_t m = P.getOptionLongValue("-m", 10*n);
  bool ordered = P.getOption("-o");
  bool adjArray = P.getOption("-j");
  edgeArray<size_t> EA = edgeRandomWithDimension<size_t>(dim, m/n, n);
  writeGraphFromEdges(EA, fname, adjArray, ordered);
}
