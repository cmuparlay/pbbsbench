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


#include "pbbslib/parse_command_line.h"
#include "common/graph.h"
#include "common/graphIO.h"

#include "common/graphUtils.h"
#include "pbbslib/parallel.h"
using namespace benchIO;
using namespace dataGen;
using namespace std;

// Star graph

template <class intT>
edgeArray<intT> star(intT m) {
  intT numCenters = 4;
  edge<intT> *E = newA(edge<intT>,m);
  parallel_for (0, m, [&] (size_t k) {
    E[k].u = utils::hash(k)%numCenters;
    E[k].v = k+numCenters;
    });
  return edgeArray<intT>(E,m+numCenters,m+numCenters,m);
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-m <numedges>] [-d <dims>] [-o] [-j] n <outFile>");
  pair<intT,char*> in = P.sizeAndFileName();
  intT n = in.first;
  char* fname = in.second;
  int dim = P.getOptionIntValue("-d", 0);
  bool ordered = P.getOption("-o");
  bool adjArray = P.getOption("-j");
  edgeArray<intT> EA = star<intT>(n);
  int r;
  if (adjArray) {
    graph<intT> G = graphFromEdges<intT>(EA,1);
    EA.del();
    if (!ordered) G = graphReorder<intT>(G);
    r = writeGraphToFile<intT>(G, fname);
    G.del();
  } else {
    graph<intT> G = graphFromEdges<intT>(EA,0);
    EA.del();
    if (!ordered) G = graphReorder<intT>(G);
    EA = edgesFromGraph<intT>(G);
    G.del();
    if (!ordered) std::random_shuffle(EA.E, EA.E + EA.nonZeros);
    r = writeEdgeArrayToFile<intT>(EA, fname);
    EA.del();
  }
  return r;
}
