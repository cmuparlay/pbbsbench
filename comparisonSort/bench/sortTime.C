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
#include "parlay/random.h"
#include "parlay/parallel.h"
#include "common/sequenceIO.h"
#include "common/parseCommandLine.h"
#include "common/time_loop.h"

using namespace std;
using namespace benchIO;

template <typename T, typename Less>
int timeSort(sequence<parlay::chars> const &In, Less less, int rounds, bool permute, char* outFile) {
  sequence<T> A = parseElements<T>(In.cut(1, In.size()));
  
  size_t n = A.size();
  if (permute) A = parlay::random_shuffle(A);
  sequence<T> B;
  time_loop(rounds, 2.0,
	    [&] () {if constexpr(INPLACE) B = A;},
	    [&] () {
	      if constexpr(INPLACE) compSort(B, less);
	      else B = compSort(A, less);},
	    [&] () {});
  cout << endl;
  if (outFile != NULL) writeSequenceToFile(B, outFile);
  return 0;
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-p] [-o <outFile>] [-r <rounds>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  int rounds = P.getOptionIntValue("-r",1);
  bool permute = P.getOption("-p");

  auto In = get_tokens(iFile);
  elementType in_type = elementTypeFromHeader(In[0]);
  size_t n = In.size() - 1;


  if (in_type == intType) {
    return timeSort<int>(In, std::less<int>(), rounds, permute, oFile);
  } else if (in_type == doubleT) {
    return timeSort<double>(In, std::less<double>(), rounds, permute, oFile);
  } else if (in_type == intPairT) {
    using ipair = pair<int,int>;
    auto less = [] (ipair a, ipair b) {return a.first < b.first;};
    return timeSort<ipair>(In, less, rounds, permute, oFile);
  } else if (in_type == doublePairT) {
    using dpair = pair<double,double>;
    auto less = [] (dpair a, dpair b) {return a.first < b.first;};
    return timeSort<dpair>(In, less, rounds, permute, oFile);
  } else if (in_type == stringT) {
    using str = parlay::chars;
    auto strless = [&] (str const &a, str const &b) -> bool {
      auto sa = a.data();
      auto sb = b.data();
      auto ea = sa + min(a.size(),b.size());
      while (sa < ea && *sa == *sb) {sa++; sb++;}
      return sa == ea ? (a.size() < b.size()) : *sa < *sb;
    };
    return timeSort<str>(In, strless, rounds, permute, oFile); 
  } else {
    cout << "sortTime: input file not of right type" << endl;
    return(1);
  }
}
