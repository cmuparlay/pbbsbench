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
#include "common/sequenceIO.h"
#include "common/atomics.h"
#include "common/parse_command_line.h"
using namespace std;
using namespace benchIO;

template <class T, class LESS>
void checkSort(sequence<sequence<char>> In,
	       sequence<sequence<char>> Out,
	       LESS less) {
  sequence<T> in_vals = parseElements<T>(In.cut(1, In.size()));
  sequence<T> out_vals = parseElements<T>(Out.cut(1, In.size()));
  size_t n = in_vals.size();
  auto sorted_in = parlay::stable_sort(in_vals, less);
  size_t error = n;
  parlay::parallel_for (0, n, [&] (size_t i) {
      if (out_vals[i] != sorted_in[i]) 
	pbbs::write_min(&error,i,std::less<size_t>());
  });
  if (error < n) {
    auto expected = parlay::to_chars(sorted_in[error]);
    auto got = parlay::to_chars(out_vals[error]);
    cout << "integer sort: check failed at location i=" << error
	 << " expected " << expected << " got " << got << endl;
    abort();
  }
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"<inFile> <outFile>");
  pair<char*,char*> fnames = P.IOFileNames();
  char* infile = fnames.first;
  char* outfile = fnames.second;
  
  auto In = get_tokens(infile);
  elementType in_type = elementTypeFromHeader(In[0]);
  size_t in_n = In.size() - 1;

  auto Out = get_tokens(outfile);
  elementType out_type = elementTypeFromHeader(Out[0]);
  size_t out_n = In.size() - 1;

  if (in_type != out_type) {
    cout << argv[0] << ": in and out types don't match" << endl;
    return(1);
  }
  
  if (in_n != out_n) {
    cout << argv[0] << ": in and out lengths don't match" << endl;
    return(1);
  }

  auto less = [&] (uint a, uint b) {return a < b;};
  auto lessp = [&] (uintPair a, uintPair b) {return a.first < b.first;};
  
  switch (in_type) {
  case intType: 
    checkSort<uint>(In, Out, less);
    break; 
  case intPairT: 
    checkSort<uintPair>(In, Out, lessp);
    break; 
  default:
    cout << argv[0] << ": input files not of right type" << endl;
    return(1);
  }
}
