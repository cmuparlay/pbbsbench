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
#include "wc.h"
using namespace std;
using namespace benchIO;

using str_t = parlay::sequence<char>;

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"<infile> <outfile>");
  pair<char*,char*> fnames = P.IOFileNames();
  str_t In = readStringFromFile(fnames.first);
  str_t Out = readStringFromFile(fnames.second);

  auto rin = wordCounts(In, false);

  auto tokens = parlay::tokens(Out, is_space);
  auto rout = parlay::tabulate(tokens.size()/2, [&] (size_t i) {
      return result_type(tokens[2*i], parlay::chars_to_long(tokens[2*i+1])); });

  if (rin.size() != rout.size()) {
    cout << "wcCheck: number of unique words do not match, " << rin.size() << " found in the input, but " << rout.size() << " found in the supplied result. " << endl;
    return(1);
  }

  auto strless = [] (str_t const &a, str_t const &b) -> bool {
    auto sa = a.data();
    auto sb = b.data();
    auto ea = sa + min(a.size(),b.size());
    while (sa < ea && *sa == *sb) {sa++; sb++;}
    return sa == ea ? (a.size() < b.size()) : *sa < *sb;
  };

  auto cmp = [&] (result_type a, result_type b) {
    return (a.second < b.second ||
	    ((a.second == b.second) && strless(a.first,b.first)));};

  rout = parlay::sort(rout,cmp);
  rin = parlay::sort(rin,cmp);
  
  if (!(rin == rout)) {
    cout << "wcCheck: counts do not match" << endl;
    return(1);
  }

  return 0;
}
