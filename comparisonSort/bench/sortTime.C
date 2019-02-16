// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2010 Guy Blelloch and the PBBS team
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
#include "get_time.h"
#include "random_shuffle.h"
#include "parallel.h"
#include "sequenceIO.h"
#include "parseCommandLine.h"

using namespace std;
using namespace benchIO;

template <class T, class CMP>
void timeSort(T* Ap, size_t n, CMP f, int rounds, bool permute, char* outFile) {
  timer t;
  sequence<T> A(Ap,n);
  if (permute) A = pbbs::random_shuffle(A);
  sequence<T> B(n);
  parallel_for (0, n, [&] (size_t i) {B[i] = A[i];});
  compSort(B.begin(), n, f); // run one sort to "warm things up"
  for (int i=0; i < rounds; i++) {
    parallel_for (0, n, [&] (size_t i) {B[i] = A[i];});
    t.start();
    compSort(B.begin(), n, f);
    t.next("");
  }
  cout << endl;
  if (outFile != NULL) writeSequenceToFile(B, outFile);
}

template <class T, class CMP>
void timeSortPair(T* Ap, size_t n, CMP f, int rounds, bool permute, char* outFile) {
  timer t;
  sequence<T> A(Ap,n);
  cout << sizeof(T) << endl;
  cout << A[0].first << endl;
  if (permute) A = pbbs::random_shuffle(A);
  sequence<T> B(n);
  parallel_for (0, n, [&] (size_t i) {B[i] = A[i];});
  compSort(B.begin(), n, f); // run one sort to "warm things up"
  for (int i=0; i < rounds; i++) {
    parallel_for (0, n, [&] (size_t i) {B[i] = A[i];});
    t.start();
    compSort(B.begin(), n, f);
    t.next("");
  }
  cout << endl;
  cout << B[0].first << endl;
  if (outFile != NULL) writeSequenceToFile(B, outFile);
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-p] [-o <outFile>] [-r <rounds>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  int rounds = P.getOptionIntValue("-r",1);
  bool permute = P.getOption("-p");
  seqData D = readSequenceFromFile(iFile);
  size_t dt = D.dt;

  using intpair = pair<long,long>;
  using doublepair = pair<double,double>;
  auto lesspair = [&] (doublePair a, doublePair b) {return a.first < b.first;};

  auto strless = [&] (char* a, char* b) {return strcmp(a,b) < 0;};
  switch (dt) {
  case intType:
    timeSort((int*) D.A, D.n, std::less<int>(), rounds, permute, oFile);
    break;
  case doubleT:
    timeSort((double*) D.A, D.n, std::less<double>(), rounds, permute, oFile);
    break;
  case intPairT:
    timeSort((intpair*) D.A, D.n, std::less<intpair>(), rounds, permute, oFile);
    break;
  case doublePairT:
    timeSort((doublepair*) D.A, D.n, lesspair, rounds, permute, oFile);
    break;    
  case stringT:
    timeSort((char**) D.A, D.n, strless, rounds, permute, oFile); 
    break;
  default:
    cout << "comparisonSort: input file not of right type" << endl;
    return(1);
  }
}
