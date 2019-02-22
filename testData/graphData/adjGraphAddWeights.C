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

// Adds a random integer weight to each edge

#include <math.h>

#include "common/graph.h"
#include "common/graphIO.h"
#include "pbbslib/parse_command_line.h"

#include "pbbslib/parallel.h"
using namespace benchIO;
using namespace dataGen;
using namespace std;

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"<inFile> <outFile>");
  pair<char*,char*> fnames = P.IOFileNames();
  char* iFile = fnames.first;
  char* oFile = fnames.second;

  graph<intT> G = readGraphFromFile<intT>(iFile);

  uintT m = G.m;
  intT n = G.n;
  
  intT* Weights = newA(intT,m);

  intT maxEdgeLen = utils::log2Up(n);
  intT* Choices = newA(intT,2*maxEdgeLen);
  
  parallel_for (0, maxEdgeLen, [&] (size_t i) {
    Choices[2*i] = i+1;
    Choices[2*i+1] = i+1;
    //Choices[2*i+1] = -(i/10)-1;
    });

  parallel_for (0, m, [&] (size_t i) {
    Weights[i] = Choices[utils::hash(i) % (2*maxEdgeLen)];
    if(i%1000==0 && Weights[i] < 0) Weights[i]*=-1;
    });
  free(Choices);
  
  wghVertex<intT>* WV = newA(wghVertex<intT>,n);
  intT* Neighbors_start = G.allocatedInplace+2+n;

  parallel_for(0, n, [&] (size_t i) {
    WV[i].Neighbors = G.V[i].Neighbors;
    WV[i].degree = G.V[i].degree;
    intT offset = G.V[i].Neighbors - Neighbors_start;
    WV[i].nghWeights = Weights+offset;
    });

  //symmetrize
  parallel_for(0, n, [&] (size_t i) {
      parallel_for(0, WV[i].degree, [&] (size_t j) {
	  intT ngh = WV[i].Neighbors[j];
	  if(ngh > i) {
	    for(intT k=0;k<WV[ngh].degree;k++) {
	      intT ngh_ngh = WV[ngh].Neighbors[k];
	      if(ngh_ngh == i) {
		WV[i].nghWeights[j] = WV[ngh].nghWeights[k];
		break;
	      }
	    }
	  }
	});
    });

  wghGraph<intT> WG(WV,n,m,(intT*)G.allocatedInplace,Weights);
  int r = writeWghGraphToFile<intT>(WG,oFile);

  WG.del();
  
  return r;
}
