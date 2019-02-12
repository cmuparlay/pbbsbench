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
#include "parallel.h"
#include "IO.h"
#include "sequenceIO.h"
#include "parse_command_line.h"
#include "sequence_ops.h"

#include "SA.h"
using namespace std;
using namespace benchIO;
using uint = unsigned int;
using uchar = unsigned char;

template<class F, class G, class H>
void loop(int rounds, F initf, G runf, H endf) {
  timer t;
  initf();
  runf();
  endf();
  for (int i=0; i < rounds; i++) {
    initf();
    t.start();
    runf();
    t.next("");
    endf();
  }
}

void timeSuffixArray(pbbs::sequence<char> const &s, int rounds, char* outFile) {
  size_t n = s.size();
  pbbs::sequence<uchar> ss(n, [&] (size_t i) {return (uchar) s[i];});
  pbbs::sequence<uint> R;
  loop(rounds,
       [&] () {R = sequence<uint>();},
       [&] () {R = suffixArray(ss);},
       [&] () {}
       );
  //cout<<"Peak-memory: " <<getPeakRSS() / (1024*1024) << endl;
  cout << endl;
  if (outFile != NULL) writeSequenceToFile(R, outFile);
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] [-r <rounds>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  int rounds = P.getOptionIntValue("-r",1);
  pbbs::sequence<char> S = readStringFromFile(iFile);
  
  timeSuffixArray(S, rounds, oFile);
}
