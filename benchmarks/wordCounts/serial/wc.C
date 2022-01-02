
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
#include <unordered_map>
#include "parlay/primitives.h"
#include "parlay/io.h"
#include "parlay/internal/get_time.h"
#include "wc.h"

using namespace std;

parlay::sequence<result_type> wordCounts(charseq const &s, bool verbose=false) {
  parlay::internal::timer t("word counts", verbose);
  if (verbose)
    cout << "number of characters = " << s.size() << endl;
  
  // copy to mutable vector and convert non alpha chars to spaces
  vector<char> str(s.size()+1);
  for (size_t i=0; i < s.size(); i++) {
    char c = s[i];
    if (c >= 65 && c < 91) str[i] = c + 32;   // upper to lower
    else if (c >= 97 && c < 123) str[i] = c;  // already lower
    else str[i] = ' ';                     // all other
  }
  str[s.size()] = 0;
  t.next("copy");
  
  // tokenize
  vector<char*> tokens;
  char* next_token = strtok(str.data(), " ");
  size_t count = 0;
  while (next_token != NULL) {
    tokens.push_back(next_token);
    next_token = strtok (NULL, "\r\n\t ");
    count++;
  }
  if (verbose)
    cout << "number of words = " << count << endl;
  t.next("tokenize");
  
  // define a hash table
  auto strhash = [] (char* s) -> size_t {
    size_t hash = 5381;
    char* l = s;
    while (*l != 0) {
      hash = ((hash << 5) + hash) + *l; l++;}
    return hash;
  };

  auto streql = [] (char* a, char* b) {
    return strcmp(a,b) == 0;};

  unordered_map<char*, size_t, decltype(strhash), decltype(streql)> word_map(count, strhash, streql);

  // add each token to the word_map
  for (size_t i=0; i < count; i++)
    ++word_map[tokens[i]];
  t.next("insert into hash table");
  if (verbose)
    cout << "result size = " << word_map.size() << endl;
  
  // pull out elements into a sequence
  parlay::sequence<result_type> result;
  result.reserve(word_map.size());
  for (const auto &pair : word_map) 
    result.push_back(result_type(parlay::to_chars(pair.first),pair.second));
  t.next("extract results");
  
  return result;
}
