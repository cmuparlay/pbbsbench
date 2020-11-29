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
#include "index.h"

namespace delayed = parlay::block_delayed;
using namespace std;

charseq build_index(charseq const &s, charseq const &doc_start) {
  parlay::internal::timer t("build Index", true);
  size_t n = s.size();
  size_t m = doc_start.size();
  //cout << "number of characters = " << n << endl;

  auto starts = delayed::filter(parlay::iota(n-m+1), [&] (size_t i) {
    for (int j=0; j < m; j++)
      if (doc_start[j] != s[i+j]) return false;
    return true;});
  auto num_docs = starts.size();
  t.next("filter index");
  cout << "num docs = " << num_docs << endl;

  auto strless = [&] (charseq const &a, charseq const &b) -> bool {
      auto sa = a.data();
      auto sb = b.data();
      auto ea = sa + min(a.size(),b.size());
      while (sa < ea && *sa == *sb) {sa++; sb++;}
      return sa == ea ? (a.size() < b.size()) : *sa < *sb;
    };

  auto x = parlay::tabulate(num_docs, [&] (size_t i) {
    size_t end = (i==num_docs-1) ? n : starts[i+1];
    auto str = parlay::map(s.cut(starts[i]+m, end), [] (char c) -> char {
	if (c >= 65 && c < 91) return c + 32;
	else if (c >= 97 && c < 123) return c;
	else return 0;});
    auto is_space = [] (char c) {return c == 0;};
    auto t = parlay::unique(parlay::sort(parlay::tokens(str, is_space), strless));
    return parlay::map(t, [&] (auto str) {
      return std::pair(str, i);});
  });
  t.next("tabulate");

  auto y = parlay::flatten(x);
  t.next("flatten");
  cout << "num words in docs = " << y.size() << endl;

  auto z = parlay::group_by_key(y, strless);
  t.next("group by");

  cout << "num unique words = " << z.size() << endl;

  charseq space(1, ' ');
  charseq newline(1, '\n');

  auto b = parlay::map(z, [&] (auto x) -> charseq {
     auto docs = x.second;
     size_t len = docs.size()*2 + 2;
     auto ss = parlay::tabulate(len, [&] (size_t i) {
       if (i == 0) return x.first;							   if (i == len-1) return newline;
       if (i%2 == 1) return space;
       return parlay::to_chars(docs[i/2-1].second);});
     return parlay::flatten(ss);});
  t.next("to strings");

  auto c = parlay::flatten(b);
  t.next("flatten");
  return c;
}
