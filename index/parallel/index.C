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
#include "parlay/io.h"
#include "parlay/internal/group_by.h"
#include "parlay/internal/get_time.h"
#include "index.h"

namespace delayed = parlay::block_delayed;
using namespace std;

charseq build_index(charseq const &s, charseq const &doc_start,
		    bool verbose = false) {
  parlay::internal::timer t("build Index", verbose);
  size_t n = s.size();
  size_t m = doc_start.size();

  // sequence of indices to the start of each document
  auto starts = delayed::filter(parlay::iota(n-m+1), [&] (size_t i) {
    for (int j=0; j < m; j++)
      if (doc_start[j] != s[i+j]) return false;
    return true;});
  auto num_docs = starts.size();
  t.next("get starts");
  if (verbose) cout << "num docs = " << num_docs << endl;

  // generate sequence of token-doc_id pairs for each document
  auto docs = parlay::tabulate(num_docs, [&] (unsigned int doc_id) {
    size_t start = starts[doc_id] + m;					
    size_t end = (doc_id==num_docs-1) ? n : starts[doc_id+1];

    // blank out all non characters, and convert to lowercase
    auto str = parlay::map(s.cut(start, end), [] (char c) -> char {
	if (c >= 65 && c < 91) return c + 32;   // upper to lower
	else if (c >= 97 && c < 123) return c;  // already lower
	else return 0;});                       // all other

    // generate tokens (i.e., contiguous regions of non-zero characters)
    auto tokens = parlay::tokens(str, [] (char c) {return c == 0;});

    // remove duplicate tokens
    tokens = parlay::remove_duplicates(std::move(tokens));

    // tag each remaining token with document id
    return parlay::map(tokens, [&] (auto str) {
        return std::pair(str, doc_id);});
  });
  t.next("generate document tokens");

  auto word_doc_pairs = parlay::flatten(std::move(docs));
  t.next("flatten document tokens");
  if (verbose)
    cout << "num words in docs = " << word_doc_pairs.size() << endl;

  // group by word, each with a sequence of docs it appears in.
  auto words = parlay::group_by_key(std::move(word_doc_pairs));
  t.next("group by word");
  if (verbose)
    cout << "num unique words = " << words.size() << endl;
  
  parlay::sort_inplace(words, [] (auto const &l, auto const &r) {
			           return l.first < r.first;});
  t.next("sort words");

  // generate string for each document number
  auto docstr = parlay::tabulate(num_docs, [] (size_t i) {
		     return parlay::to_chars(i);});
  
  // format output for each word
  charseq space(1, ' ');
  charseq newline(1, '\n');
  auto b = parlay::map(words, [&] (auto const &wd_pair) -> charseq {
     //auto [word, doc_ids] = std::move(wd_pair);
     auto word = std::move(wd_pair.first);
     auto doc_ids = std::move(wd_pair.second);
     size_t len = doc_ids.size()*2 + 2;
     // each line consists of the word followed by
     // the list of documents ids separared by spaces 
     // and terminated by a newline.
     auto ss = parlay::tabulate(len, [&] (size_t i) {
       if (i == 0) return word;
       if (i == len-1) return newline;
       if (i%2 == 1) return space;
       return docstr[doc_ids[i/2-1]];});
     return parlay::flatten(ss);});
  t.next("format words");

  // flatten across words
  auto c = parlay::flatten(std::move(b));
  t.next("flatten formatted words");
  return c;
}
