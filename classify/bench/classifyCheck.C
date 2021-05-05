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
using namespace std;
using namespace benchIO;

using str_t = parlay::sequence<char>;
using value = unsigned char;
using row = parlay::sequence<value>;

auto read_row(string filename) {
  auto is_item = [] (char c) -> bool { return c == ',';};
  auto str = parlay::chars_from_file(filename);
  return parlay::map(parlay::tokens(str, is_item), 
		     [] (auto s) -> value {return parlay::chars_to_int(s);});
}

void report_correct(row result, row labels) {
  size_t n = result.size();
  if (n != labels.size()) {
    cout << "size mismatch of results and labels" << endl;
    return;
  }
  size_t num_correct = parlay::reduce(parlay::tabulate(n, [&] (size_t i) {
         return (result[i] == labels[i]) ? 1 : 0;}));
  float percent_correct = (100.0 * num_correct)/n;
  cout << num_correct << " correct out of " << n
       << ", " << percent_correct << " percent" << endl;
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"<infile> <outfile>");
  pair<char*,char*> fnames = P.IOFileNames();
  bool verbose = P.getOption("-v");
  string expected_file_prefix = fnames.first;
  string given_file = fnames.second;
  auto expected_labels = read_row(expected_file_prefix.append(".labels"));
  auto given_labels = read_row(given_file);
  if (verbose) report_correct(expected_labels, given_labels);

  return 0;
}
