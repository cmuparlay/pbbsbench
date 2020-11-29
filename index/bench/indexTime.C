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
#include "parlay/internal/get_time.h"
#include "common/IO.h"
#include "common/sequenceIO.h"
#include "common/parse_command_line.h"

#include "index.h"
using namespace std;
using namespace benchIO;

template<class F, class G, class H>
void loop(int rounds, bool cold, F initf, G runf, H endf) {
  parlay::internal::timer t;
  double start_time = t.get_time();
  if (!cold)
    do { // run for a couple seconds to "warm things up"
      initf(); runf(); endf();
    } while (t.get_time() < start_time + 1.0);
  for (int i=0; i < rounds; i++) {
    initf();
    t.start();
    runf();
    t.next("");
    endf();
  }
}

void timeWordCounts(charseq const &s, charseq const &start, 
		    int rounds, bool cold, char* outFile) {
  size_t n = s.size();
  charseq R;
  loop(rounds, cold,
       [&] () {R.clear();},
       [&] () {R = build_index(s, start);},
       [&] () {}
       );
  cout << endl;
  if (outFile != NULL) parlay::chars_to_file(R, outFile);
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-c] [-o <outFile>] [-r <rounds>] <inFile>");
  char* iFile = P.getArgument(0);
  bool cold = P.getOption("-c");
  char* oFile = P.getOptionValue("-o");
  int rounds = P.getOptionIntValue("-r",1);
  //parlay::sequence<char> S = parlay::chars_from_file(iFile, true);
  parlay::sequence<char> S = parlay::to_sequence(parlay::file_map(iFile));
  
  string header = "<doc";
  timeWordCounts(S, parlay::to_sequence(header), rounds, cold, oFile);
}
