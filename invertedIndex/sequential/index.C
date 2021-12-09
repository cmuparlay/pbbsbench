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

#include <map>
#include <unordered_set>
#include <unordered_map>
#include <string>

#include "parlay/io.h"
#include "parlay/sequence.h"

#include "parlay/internal/get_time.h"

#include "index.h"

namespace delayed = parlay::block_delayed;
using namespace std;

charseq build_index(charseq const &s, charseq const &doc_start,
		    bool verbose = false) {
  parlay::internal::timer t("build Index", verbose);
  size_t n = s.size();
  size_t m = doc_start.size();
  
  // group by word, each with a sequence of docs it appears in.
  std::unordered_map<charseq, std::vector<size_t>> words;
  std::vector<charseq> doc_id_str;

  // Find the first document delimiter
  auto doc_begin = std::search(std::begin(s), std::end(s), std::begin(doc_start), std::end(doc_start));

  // Generate, for each document, the tokens contained within it
  for (size_t doc_id = 0; doc_begin != std::end(s); ++doc_id) {
    doc_id_str.push_back(parlay::to_chars(doc_id));
    
    // Find the end of the current document
    auto doc_end = std::search(doc_begin + m, std::end(s), std::begin(doc_start), std::end(doc_start));
    
    // generate tokens (i.e., contiguous regions of alphabetic characters)
    auto token_begin = doc_begin + m;
    while (token_begin != doc_end && !std::isalpha(*token_begin)) token_begin++;
    while (token_begin != doc_end) {
      auto token_end = std::find_if(token_begin, doc_end, [](unsigned c) { return !std::isalpha(c); });
      
      size_t token_length = token_end - token_begin;
      charseq token(token_length);
      std::transform(token_begin, token_end, std::begin(token), [](unsigned c) { return std::tolower(c); });

      auto& doc_ids_for_token = words[std::move(token)];
      if (doc_ids_for_token.empty() || doc_ids_for_token.back() != doc_id) {
        doc_ids_for_token.push_back(doc_id);
      }
      
      token_begin = token_end;
      while (token_begin != doc_end && !std::isalpha(*token_begin)) token_begin++;
    }
    
    doc_begin = doc_end;
  }
  t.next("generate document tokens");

  if (verbose)
    cout << "num unique words = " << words.size() << endl;

  std::vector<charseq> sorted_words;
  sorted_words.reserve(words.size());
  for (const auto& [word, doc_ids] : words) {
    sorted_words.push_back(word);
  }
  std::sort(std::begin(sorted_words), std::end(sorted_words));

  charseq result;
  for (const auto& word: sorted_words) {
    const auto& doc_ids = words[word];
    result.append(std::begin(word), std::end(word));
    result.push_back(' ');
    for (size_t i = 0; i < doc_ids.size(); i++) {
      result.append(doc_id_str[doc_ids[i]]);
      if (i < doc_ids.size() - 1) result.push_back(' ');
    }
    result.push_back('\n');
  }
  t.next("format words");

  return result;
}
