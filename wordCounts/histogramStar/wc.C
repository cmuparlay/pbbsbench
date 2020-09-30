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
#include "parlay/parallel_io.h"
#include "parlay/internal/collect_reduce.h"
#include "common/get_time.h"
#include "wc.h"

using namespace std;

parlay::sequence<result_type> wordCounts(charseq const &s) {
  timer t("word counts");
  cout << "number of characters = " << s.size() << endl;
  
  auto is_space = [] (char c) {
    switch (c)  {
    case '\r': 
    case '\t': 
    case '\n': 
    case 0:
    case ' ' : return true;
    default : return false;
    }
  };

  size_t n = s.size();
  charseq str = parlay::sequence<char>::from_function(n + 1, [&] (size_t i) -> char {
      return (i == n || is_space(s[i])) ? 0 : s[i];});
  t.next("copy");
  
  auto get_word = [&] (parlay::slice<char*,char*> x) -> char* {
    return x.begin();
  };

  auto words = parlay::tokens(str, [] (char c) {return c == 0;}, get_word);
  
  t.next("tokens");
  cout << "number of words = " << words.size() << endl;

  struct hasheq {
    // a simple hash function on char sequences
    static inline size_t hash(char* a) {
      size_t hash = 5381;
      for (size_t i = 0; a[i] != 0; i++) 
  	hash = ((hash << 5) + hash) + a[i];
      return hash;
    }
    static inline bool eql(char* a, char* b) {
      return strcmp(a,b) == 0; }
  };
  
  auto hist = parlay::internal::histogram_sparse(make_slice(words), hasheq());
  t.next("collect reduce");

  cout << "result.size(): " << hist.size() << endl;
  cout << "result[0]: " << hist[0].first << ", " << hist[0].second << endl;
  parlay::sequence<result_type> result =
    parlay::tabulate(hist.size(), [&] (size_t i) {
	size_t len = strlen(hist[i].first);
	return result_type(parlay::sequence<char>(std::move(hist[i].first), hist[i].first+len),
			   hist[i].second);
      });
  hist.clear_uninitialized();
  return result;
}
