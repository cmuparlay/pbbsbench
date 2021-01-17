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

#include "classify.h"

using namespace std;
using namespace benchIO;

void timeClassify(features const &Train, rows const &Test, int rounds, bool verbose, char* outFile) {
  time_loop(rounds, 2.0,
	    [&] () {},
	    [&] () {classify(Train, Test, verbose);},
	    [&] () {});
  cout << endl;
  //if (outFile != NULL) writeHistogramsToFile(R, outFile);
}

auto read_data(string filename) {
  int num_buckets = 64;
  auto is_line = [] (char c) -> bool { return c == '\r' || c == '\n'|| c == 0;};
  auto is_item = [] (char c) -> bool { return c == ',';};
  auto str = parlay::file_map(filename);
  size_t j = find_if(str,is_line) - str.begin();
  auto head = parlay::make_slice(str).cut(0,j);
  auto rest = parlay::make_slice(str).cut(j+1, str.end()-str.begin());
  auto types = parlay::map(parlay::tokens(head, is_item), [] (auto str) {return str[0];});

  auto process_line = [&] (auto line) {
      auto to_int = [&] (auto x) -> value {
	  int v = parlay::chars_to_int(parlay::to_sequence(x));
	  if (v < 0 || v > max_value) {
	    cout << "entry out range: value = " << v << endl;
	    return 0;
	  }
	  return v;};
      return parlay::map_tokens(line, to_int, is_item);};
  
  return std::pair(types, parlay::map_tokens(rest, process_line, is_line));
}

auto rows_to_features(sequence<char> types, rows const &A) {
  int num_features = types.size();
  size_t num_rows = A.size();

  auto get_feature = [&] (size_t io) {
    int i = (io == 0) ? num_features - 1 : io - 1; // put label at front
    bool is_discrete = (types[i] == 'd');
    auto vals = parlay::tabulate(num_rows, [&] (size_t j) -> value {return A[j][i];});
    int maxv = parlay::reduce(vals, parlay::maxm<char>());
    return feature(is_discrete, maxv+1, vals);
  };

  return parlay::tabulate(num_features, get_feature);
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] [-r <rounds>] <inFile>");
  string iFile = P.getArgument(0);
  string train_file = iFile;
  string test_file = iFile;
  char* oFile = P.getOptionValue("-o");
  bool verbose = P.getOption("-v");
  int rounds = P.getOptionIntValue("-r",1);
  auto [types, Train] = read_data(train_file.append(".train"));
  auto [xtypes, Test] = read_data(test_file.append(".test"));
  cout << "num tests: " << Test.size() << endl;
  features TrainFeatures = rows_to_features(types, Train);
  timeClassify(TrainFeatures, Test, rounds, verbose, oFile);
}
