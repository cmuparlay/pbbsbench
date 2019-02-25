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

#include "pbbslib/parallel.h"
#include "pbbslib/parse_command_line.h"
#include "pbbslib/get_time.h"
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
template <class intT>
edgeArray<intT> edgeRandomWithDimension(intT dim, intT degree, intT numRows) {
  intT nonZeros = numRows*degree;
  pbbs::sequence<edge<intT>> E(nonZeros, [&] (size_t k) {
      intT i = k / degree;
      intT j;
      if (dim==0) {
	intT h = k;
	do {
	  j = ((h = dataGen::hash<intT>(h)) % numRows);
	} while (j == i);
      } else {
	intT pow = dim+2;
	intT h = k;
	do {
	  while ((((h = dataGen::hash<intT>(h)) % 1000003) < 500001)) pow += dim;
	  j = (i + ((h = dataGen::hash<intT>(h)) % (((long) 1) << pow))) % numRows;
	} while (j == i);
      }
      return edge<intT>(i, j);
    });
  return edgeArray<intT>(std::move(E), numRows, numRows);
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-m <numedges>] [-d <dims>] [-o] [-j] n <outFile>");
  pair<intT,char*> in = P.sizeAndFileName();
  intT n = in.first;
  char* fname = in.second;
  int dim = P.getOptionIntValue("-d", 0);
  intT m = P.getOptionLongValue("-m", 10*n);
  bool ordered = P.getOption("-o");
  bool adjArray = P.getOption("-j");
  edgeArray<intT> EA = edgeRandomWithDimension<intT>(dim, m/n, n);
  writeGraphFromEdges(EA, fname, adjArray, ordered);
}
