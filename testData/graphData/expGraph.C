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
#include "sequence.h"
using namespace benchIO;
using namespace dataGen;
using namespace std;

template <class intT>
edgeArray<intT> expDegrees(intT n) {
  intT* Degrees = newA(intT,n);
  parallel_for(intT i=0;i<n;i++){
    intT rand = dataGen::hash<intT>(i)%n;
    intT p = 1;
    while(1) {
      if(rand < p) {
	Degrees[i] = max(1,n/p);
	break;
      }
      p*=2;
    } 
  }

  intT totalDegree = sequence::plusScan(Degrees,Degrees,n);
  edge<intT> *E = newA(edge<intT>,totalDegree);
  parallel_for(intT ii=0;ii<n;ii++){
    intT start = Degrees[ii];
    intT end = (ii < n-1) ? Degrees[ii+1] : totalDegree;
    if(end-start > 1000)
      parallel_for(intT j=start;j<end;j++){
	E[j].u = ii; E[j].v = utils::hash2(j) % n;
      }
    else
      for(intT j=start;j<end;j++){
	E[j].u = ii; E[j].v = utils::hash2(j) % n;
      }
  }
  free(Degrees);
  return edgeArray<intT>(E,n,n,totalDegree);
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
  edgeArray<intT> EA = expDegrees<intT>(n);
  int r;
  if (adjArray) {
    graph<intT> G = graphFromEdges<intT>(EA,1);
    EA.del();
    if (!ordered) G = graphReorder<intT>(G);
    r = writeGraphToFile<intT>(G, fname);
    G.del();
  } else {
    if (!ordered) std::random_shuffle(EA.E, EA.E + EA.nonZeros);
    r = writeEdgeArrayToFile<intT>(EA, fname);
    EA.del();
  }
  return r;
}
