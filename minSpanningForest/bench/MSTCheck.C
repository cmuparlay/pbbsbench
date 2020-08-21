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
#include <cstring>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/IO.h"
#include "common/graph.h"
#include "common/graphIO.h"
#include "common/parse_command_line.h"
#include "MST.h"

using namespace std;


int main(int argc, char* argv[]) { 
  commandLine P(argc,argv,"<inFile> <outfile>");
  pair<char*,char*> fnames = P.IOFileNames();
  char* iFile = fnames.first;
  char* oFile = fnames.second;
  using WE = wghEdge<vertexId,edgeWeight>;
  using graph = wghEdgeArray<vertexId,edgeWeight>;
  
  graph In = readWghEdgeArrayFromFile<vertexId,edgeWeight>(iFile);
  parlay::sequence<edgeId> Out = readIntSeqFromFile<edgeId>(oFile);
  size_t n = Out.size();
  size_t in_m = In.m;
  //check num edges
  parlay::sequence<edgeId> serialMST = mst(In);
  if (n != serialMST.size()) {
    cout << "Wrong edge count: MST has " << serialMST.size()
	 << " edges but algorithm returned " << n << " edges\n";
    return (1);
    }
  
  //check for cycles
  parlay::sequence<bool> flags(in_m, false);
  parlay::parallel_for(0, n, [&] (size_t i) {flags[Out[i]] = true;});
  
  parlay::sequence<WE> E = parlay::pack(In.E, flags);
  graph EA(std::move(E), In.n);
  
  parlay::sequence<edgeId> check = mst(EA);
  if (n != check.size()){
    cout << "Result is not a spanning forest : it has a cycle" << endl;
    return (1);
  }
  
  // check total weight
  double total1 = parlay::reduce(parlay::delayed_seq<double>(n, [&] (size_t i) {
	return In[Out[i]].weight;}), parlay::addm<double>());

  double total2 = parlay::reduce(parlay::delayed_seq<double>(n, [&] (size_t i) {
	return In[serialMST[i]].weight;}), parlay::addm<double>());

  // need to account for roundoff error, and possible equal weight trees
  if (fabs((total1 - total2)/total1) > 1e-9) {
    cout << "MST has a weight of " << total1 <<
      " but should have a weight of " << total2 << endl;
    return (1);
  }

  return 0;
}
