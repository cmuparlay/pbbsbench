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
#include <cstring>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/IO.h"
#include "common/parse_command_line.h"
#include "common/atomics.h"

using namespace std;
using namespace benchIO;

typedef unsigned char uchar;

bool strLessBounded (const uchar* s1, const uchar* s2, long n, const uchar* end) {
  while (*s1==*s2) {
    if (n-- < 0) return 1;
    if (++s1 == end) return 1;
    if (++s2 == end) return 0;
  };
  return (*s1 < *s2);
}

bool isPermutation(parlay::sequence<long> const &SA) {
  size_t n = SA.size();
  parlay::sequence<long> seen(n,(long) 0);
  parlay::parallel_for (0, n, [&] (size_t i) {seen[SA[i]] = 1;});
  long nseen = parlay::reduce(seen, parlay::addm<long>());
  return (nseen == n);
}

bool isSorted(parlay::sequence<uchar> const &s, parlay::sequence<long> const &SA) {
  auto p = s.begin();
  size_t n = s.size();
  int checkLen = 100;
  size_t error = n;
  parlay::parallel_for (0, n-1, [&] (size_t i) {
      if (!strLessBounded(p+SA[i], p+SA[i+1], checkLen, p + n)) {
	pbbs::write_min(&error,i,std::less<size_t>());
      }
    });
  if (error != n) {
    cout << "Suffix Array Check: not sorted at i = " 
	 << error+1 << endl;
    return 0;
  }
  return 1;
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"<infile> <outfile>");
  pair<char*,char*> fnames = P.IOFileNames();
  parlay::sequence<char> InX = readStringFromFile(fnames.first);
  auto In = parlay::tabulate(InX.size(), [&] (size_t i) -> uchar {
      return (uchar) InX[i];});
  parlay::sequence<long> Out = readIntSeqFromFile<long>(fnames.second);
  if (In.size() != Out.size()) {
    cout << "Suffix Array Check: lengths don't match" << endl;
    return 1;
  }
  if (!isPermutation(Out)) {
    cout << "Suffix Array Check: array is not a permutation" << endl;
    return 1;
  }
  return !isSorted(In, Out);
}
