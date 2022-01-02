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
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/io.h"
#include "parlay/internal/group_by.h"
#include "parlay/internal/get_time.h"
#include "wc.h"

using namespace std;

parlay::sequence<result_type> wordCounts(charseq const &s, bool verbose=false) {
  parlay::internal::timer t("word counts", verbose);
  if (verbose) cout << "number of characters = " << s.size() << endl;

  // blank out all non alpha characters, and convert upper to lowercase
  auto str = parlay::map(s, [] (char c) -> char {
    if (c >= 65 && c < 91) return c + 32;   // upper to lower
    else if (c >= 97 && c < 123) return c;  // already lower
    else return 0;});                       // all other
  t.next("copy");
  
  // generate tokens (i.e., contiguous regions of non-zero characters)
  auto words = parlay::tokens(str, [] (char c) {return c == 0;});
  t.next("tokens");
  if (verbose) cout << "number of words = " << words.size() << endl;

  auto result = parlay::histogram_by_key(std::move(words)); 
  t.next("count by key");
  cout << words.size() << endl;

  words.clear();
  t.next("clear");

  if (verbose) cout << "distinct words: " << result.size() << endl;
  return result;
}
