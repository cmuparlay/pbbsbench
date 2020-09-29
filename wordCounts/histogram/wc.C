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
  auto is_space = [] (char c) {return std::isspace(c);};
  auto words = parlay::tokens(s, is_space);
  t.next("tokens");
  cout << "number of words = " << words.size() << endl;

  struct hasheq {
    // a simple hash function on char sequences
    static inline size_t hash(charseq const &a) {
      size_t hash = 5381;
      for (size_t i = 0; i < a.size(); i++) 
	hash = ((hash << 5) + hash) + a[i];
      return hash;
    }
    static inline bool eql(charseq const &a, charseq const &b) {
      return a == b; }
  };
  
  auto result = parlay::internal::histogram_sparse(make_slice(words), hasheq());
  t.next("collect reduce");

  cout << "result.size(): " << result.size() << endl;
  cout << "result[0]: " << result[0].first << ", " << result[0].second << endl;
  return result;
}
