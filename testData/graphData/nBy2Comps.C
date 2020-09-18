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
#include "parlay/primitives.h"
#include "common/parse_command_line.h"
#include "common/graph.h"
#include "common/graphIO.h"
#include "common/graphUtils.h"
using namespace benchIO;
using namespace dataGen;
using namespace std;

template <class intV>
edgeArray<intV> genGraph(size_t n) {
  size_t m = n/2;
  auto E = parlay::tabulate(m, [&] (size_t i) -> edge<intV> {
      return edge<intV>(2*i, 2*i + 1); });
  return edgeArray<intV>(E,n,n);
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-j] [-o] n <outFile>");
  pair<int,char*> in = P.sizeAndFileName();
  size_t n = in.first;
  char* fname = in.second;
  bool ordered = P.getOption("-o");
  bool adjArray = P.getOption("-j");
  if (n % 2 != 0) {
    cout << "n must be even" << endl;
    P.badArgument();
  }
  auto EA = genGraph<size_t>(n);
  writeGraphFromEdges(EA, fname, adjArray, ordered);
}
