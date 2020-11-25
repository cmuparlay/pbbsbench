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

// SA.h defines indexT, which is the type of integer used for the elements of the
// suffix array
#include "wc.h"
using namespace std;
using namespace benchIO;

void writeHistogramsToFile(parlay::sequence<result_type> const results, char* outFile) {
  auto space = parlay::to_chars(' ');
  auto newline = parlay::to_chars('\n');
  auto str = parlay::flatten(parlay::map(results, [&] (result_type x) {
	sequence<sequence<char>> s = {
	  x.first, space, parlay::to_chars(x.second), newline};
	return flatten(s);}));
  parlay::chars_to_file(str, outFile);
}
  
template<class F, class G, class H>
void loop(int rounds, F initf, G runf, H endf) {
  parlay::internal::timer t;
  double start_time = t.get_time();
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

void timeWordCounts(parlay::sequence<char> const &s, int rounds, char* outFile) {
  size_t n = s.size();
  parlay::sequence<result_type> R;
  loop(rounds,
       [&] () {R.clear();},
       [&] () {R = wordCounts(s);},
       [&] () {}
       );
  cout << endl;
  if (outFile != NULL) writeHistogramsToFile(R, outFile);
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] [-r <rounds>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  int rounds = P.getOptionIntValue("-r",1);
  parlay::sequence<char> S = parlay::chars_from_file(iFile, true);
  //parlay::sequence<char> S = parlay::to_sequence(parlay::file_map(iFile));
  
  timeWordCounts(S, rounds, oFile);
}
