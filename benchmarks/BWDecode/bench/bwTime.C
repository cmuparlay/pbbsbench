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
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/io.h"
#include "common/time_loop.h"
#include "common/IO.h"
#include "common/sequenceIO.h"
#include "common/parse_command_line.h"
#include "algorithm/bw_encode.h"

#include "bw.h"
using namespace std;
using namespace benchIO;

auto timeBW(ucharseq const &s, int rounds, char* outFile) {
  size_t n = s.size();
  ucharseq R;
  time_loop(rounds, 1.0,
       [&] () {R.clear();},
       [&] () {R = bw_decode(s);},
       [&] () {});
  cout << endl;
  if (outFile != NULL) 
    parlay::chars_to_file(parlay::map(R, [] (uchar x) {return (char) x;}), outFile);
  return R;
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] [-r <rounds>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  int rounds = P.getOptionIntValue("-r",1);
  auto S = parlay::file_map(iFile);
  //parlay::sequence<char> S = parlay::chars_from_file(iFile, true);
  auto ss = parlay::map(S, [] (char x) {return (uchar) x;});
  auto bwseq = bw_encode<unsigned int>(ss);
  
  auto R = timeBW(bwseq, rounds, oFile);
  if (R != ss) {
    cout << "bad output for bw decode" << endl;
    abort();
  }
}
