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
#include "parlay/internal/collect_reduce.h"
#include "common/get_time.h"
#include "wc.h"

using namespace std;

// this is an "optimized" version that uses pointers to the tokens rather than
// the strings themselves.   It does not seem to be any faster (or perhaps by just
// a couple percent).
parlay::sequence<result_type> wordCounts(charseq const &s, bool verbose=false) {
  timer t("word counts", verbose);
  size_t n = s.size();
  cout << "number of characters = " << n << endl;

  // pad with extra null so last string is null terminated
  auto str = parlay::tabulate(n+1, [&] (size_t i) -> char {
      char c = (i==n) ? 0 : s[i];
      if (c >= 65 && c < 91) return c + 32;   // upper to lower
      else if (c >= 97 && c < 123) return c;  // already lower
      else return 0;});                       // all other
  t.next("copy");
  
  auto get_word = [&] (parlay::slice<char*, char*> x) -> char* {
    return x.begin();
  };

  auto words = parlay::map_tokens(str, get_word, [] (char c) {return c == 0;});
  
  t.next("tokens");
  cout << "number of words = " << words.size() << endl;

  // a simple hash function on char sequences
  auto strhash = [] (char* a) {
    size_t hash = 5381;
    for (size_t i = 0; a[i] != 0; i++) 
      hash = ((hash << 5) + hash) + a[i];
    return hash;
  };
    
  auto eql = [] (char* a, char* b) {return strcmp(a,b) == 0;};

  auto hist = parlay::histogram_by_key(words, strhash, eql);
  t.next("count by key");

  cout << "distinct words: " << hist.size() << endl;
  auto result = parlay::tabulate(hist.size(), [&] (size_t i) {
      size_t len = strlen(hist[i].first);
      auto [start, count] = hist[i];
      return result_type(parlay::sequence<char>(start, start + len),
			 count);
    });

  t.next("format out");
  return result;
}
