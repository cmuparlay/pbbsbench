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
#include "parallel.h"
#include "merge_sort.h"
#include "sequenceIO.h"
#include "parse_command_line.h"
using namespace std;
using namespace benchIO;

template <class T, class LESS>
void checkIntegerSort(void* In, void* Out, size_t n, LESS less) {
  T* A = (T*) In;
  T* B = (T*) Out;
  T* Tmp = pbbs::new_array_no_init<T>(n,1);
  pbbs::merge_sort(sequence<T>(A,n), sequence<T>(Tmp,n), less, 1);
  long j = -1;
  parallel_for (0, n, [&] (size_t i) {
      if (A[i] != B[i]) j = i;});
  if (j >= 0) {
    cout << "integer sort: check failed at i=" << j << endl;
    abort();
  }
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"<inFile> <outFile>");
  pair<char*,char*> fnames = P.IOFileNames();
  seqData In = readSequenceFromFile(fnames.first);
  seqData Out = readSequenceFromFile(fnames.second);
  size_t n = In.n;
  elementType dt = In.dt;
  if (dt != Out.dt) {
    cout << argv[0] << ": types don't match" << endl;
    return(1);
  }
  if (n != Out.n) {
    cout << argv[0] << ": lengths dont' match" << endl;
    return(1);
  }
  using uint = unsigned int;
  using upair = pair<uint,int>;
  auto less = [&] (uint a, uint b) {return a < b;};
  auto lessp = [&] (upair a, upair b) {return a.first < b.first;};
  
  switch (dt) {
  case intType:
    checkIntegerSort<uint>(In.A, Out.A, n, less); break;
  case intPairT:
    checkIntegerSort<upair>(In.A, Out.A, n, lessp); break;
  default:
    cout << argv[0] << ": input files not of right type" << endl;
    return(1);
  }
}
