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

#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/time_loop.h"
#include "common/parse_command_line.h"
#include "common/sequenceIO.h"
#include <iostream>
#include <algorithm>
using namespace std;
using namespace benchIO;

template <class T>
void timeIntegerSort(sequence<sequence<char>> In, int rounds, int bits, char* outFile) {
  auto in_vals = parseElements<T>(In.cut(1, In.size()));
  size_t n = in_vals.size();
  sequence<T> R;
  time_loop(rounds, 1.0,
       [&] () {R.clear();},
       [&] () {R = int_sort(make_slice(in_vals.data(),in_vals.data()+n), bits);},
       [] () {});
  if (outFile != NULL) writeSequenceToFile(R, outFile);
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] [-r <rounds>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  int rounds = P.getOptionIntValue("-r",1);
  int bits = P.getOptionIntValue("-b",0);

  auto In = get_tokens(iFile);
  elementType in_type = elementTypeFromHeader(In[0]);
  cout << "bits = " << bits << endl;

  switch (in_type) {
  case intType: 
    timeIntegerSort<uint>(In, rounds, bits, oFile);
    break;
  case intPairT: 
    timeIntegerSort<uintPair>(In, rounds, bits, oFile);
    break;
  default:
    cout << "integer Sort: input file not of right type" << endl;
    return(1);
  }
}

