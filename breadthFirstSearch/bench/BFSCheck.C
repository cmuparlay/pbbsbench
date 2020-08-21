// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011-2019 Guy Blelloch and the PBBS team
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
#include "common/IO.h"
#include "common/graph.h"
#include "common/graphIO.h"
#include "common/parse_command_line.h"
#include "BFS.h"
using namespace std;
using namespace benchIO;

using longseq = parlay::sequence<long>;

size_t levelNumber(size_t start, size_t level, longseq &P, longseq &L, Graph &T) {
  if (L[start] != -1) {
    cout << "BFSCheck: not a tree" << endl;
    return 1;
  }
  L[start] = level;
  for (size_t i=0; i < T[start].degree; i++) {
    size_t j = T[start].Neighbors[i];
    P[j] = start;
    if (levelNumber(j, level+1, P, L, T)) return 1;
  }
  return 0;
}

// Checks if T is valid BFS tree relative to G starting at i
int checkBFS(size_t start, Graph &G, Graph &T) {
  size_t n = G.n;
  if (n != T.n) {
    cout << "BFSCheck: vertex counts don't match: " << G.n << ", " << T.n << endl;
    return 1;
  }
  if (T.m > G.n - 1) { 
     cout << "BFSCheck: too many edges in tree " << endl;
     return 1;
  }
  parlay::sequence<long> P(n, (long) -1);
  parlay::sequence<long> L(n, (long) -1);

  if (levelNumber(start, 0, P, L, T)) return 1;
  int error = 0;
  parlay::parallel_for (0, G.n, [&] (size_t i) {
      bool Check=0;
      if (L[i] == -1) {
	for (size_t j=0; j < G[i].degree; j++) {
	  size_t ngh = G[i].Neighbors[j];
	  if (L[ngh] != -1) error = 1;
	}
      } else {
	for (size_t j=0; j < G[i].degree; j++) {
	  size_t ngh = G[i].Neighbors[j];
	  if (P[i] == ngh) Check = 1;
	  else if (L[ngh] > L[i] + 1 || L[ngh] < L[i] - 1) error = 2;
	}
	if (i != start && Check == 0) error = 3;
      }
    });
  if (error == 1) cout << "BFSCheck: connected vertex not in tree " << endl;
  else if (error == 2) cout << "BFSCheck: edge spans two levels " << endl;
  else if (error == 3) cout << "BFSCheck: parent not an edge " << endl;
  else return 0;

  return 1;
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"<inFile> <outfile>");
  pair<char*,char*> fnames = P.IOFileNames();
  char* iFile = fnames.first;
  char* oFile = fnames.second;

  Graph G = readGraphFromFile<vertexId,edgeId>(iFile);
  Graph T = readGraphFromFile<vertexId,edgeId>(oFile);

  return checkBFS(0, G, T);
}
