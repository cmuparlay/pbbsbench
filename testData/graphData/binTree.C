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

#include "IO.h"
#include "parseCommandLine.h"
#include "graph.h"
#include "graphIO.h"
#include "dataGen.h"
#include "graphUtils.h"
#include "parallel.h"
#include "sequence.h"
using namespace benchIO;
using namespace dataGen;
using namespace std;

template <class intT>
edgeArray<intT> binTree(intT n) {
  intT m = n-1;
  edge<intT> *E = newA(edge<intT>,m);
  parallel_for(intT i=0;i<n/2;i++){
    intT id = i+1;
    E[2*i].u = id-1;
    E[2*i].v = 2*id-1;
    E[2*i+1].u = id-1;
    E[2*i+1].v = 2*id;
  }
  return edgeArray<intT>(E,n,n,m);
}

int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-m <numedges>] [-d <dims>] [-o] [-j] n <outFile>");
  pair<intT,char*> in = P.sizeAndFileName();
  intT n = in.first;
  char* fname = in.second;
  int dim = P.getOptionIntValue("-d", 0);
  intT m = P.getOptionLongValue("-m", 10*n);
  bool ordered = P.getOption("-o");
  bool adjArray = P.getOption("-j");
  edgeArray<intT> EA = binTree<intT>((1<<n)-1);
  int r;
  if (adjArray) {
    graph<intT> G = graphFromEdges<intT>(EA,1);
    EA.del();
    if (!ordered) G = graphReorder<intT>(G, NULL);
    r = writeGraphToFile<intT>(G, fname);
    G.del();
  } else {
    if (!ordered) std::random_shuffle(EA.E, EA.E + EA.nonZeros);
    r = writeEdgeArrayToFile<intT>(EA, fname);
    EA.del();
  }
  return r;
}
