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
#include "common/graph.h"
#include "common/graphIO.h"
#include "common/parse_command_line.h"
#include "ST.h"
using namespace std;
using namespace benchIO;

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"<inFile> <outfile>");
  pair<char*,char*> fnames = P.IOFileNames();
  char* iFile = fnames.first;
  char* oFile = fnames.second;

  edgeArray<vertexId> In = readEdgeArrayFromFile<vertexId>(iFile);
  parlay::sequence<edgeId> Out = readIntSeqFromFile<edgeId>(oFile);
  size_t n = Out.size();

  //run serial ST
  parlay::sequence<edgeId> serialST = st(In);
  if (n != serialST.size()){
    cout << "Wrong edge count: ST has " << serialST.size()
	 << " edges but algorithm returned " << n << " edges\n";
    return (1);
  }
  
  //check if ST has cycles by running serial ST on it
  //and seeing if result changes
  parlay::sequence<bool> flags(In.nonZeros, false);
  parlay::parallel_for(0, n, [&] (size_t i) {
      flags[Out[i]] = true;});
  parlay::sequence<edge<vertexId>> E = parlay::pack(In.E, flags);
  size_t m = E.size();
  
  edgeArray<vertexId> EA(std::move(E), In.numRows, In.numCols); 
  parlay::sequence<edgeId> check = st(EA);

  if (m != check.size()){
    cout << "Result is not a spanning tree " << endl;
    return (1);
  }

  return 0;
}
