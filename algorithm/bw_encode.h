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

// This code generates the burrows wheeler transform for a string s.
// It first pads the string with a null character at the front.
// It is backwards from how normally described: i.e. it sorts all suffixes, 
// and then outputs the character before the suffix in the sorted order, 
// along with the location of the first character

// example 
//  abcabcab
//  $abcabcab (after padding)
// sorted suffixes with previous character
//   b $abcabcab (using last char on wraparound)
//   c ab
//   c abcab
//   $ abcabcab  
//   a b
//   a bcab
//   a bcabcab
//   b cab
//   b cabcab
// hence the output string is: "bcc$aaabb"
#include <iostream>
#include "../parlay/parallel.h"
#include "../parlay/primitives.h"
#include "../parlay/io.h"
#include "suffix_array.h"

using uchar = unsigned char;
using ucharseq = parlay::sequence<uchar>;

// Int needs to be big enough to represent the lenght of s
// Can be unsigned
template <class Int>
ucharseq bw_encode(ucharseq const &s) {
  size_t n = s.size();

  // pad with a null at the start
  auto ss = parlay::tabulate(n+1, [&] (size_t i) {
      return i == 0 ? 0 : s[i-1];});

  // Sort on suffixes
  // for the example: <0, 7, 4, 1, 8, 5, 2, 6, 32>
  // zero will always be at the front
  auto sa = suffix_array<Int>(ss);
  // std::cout << parlay::to_chars(sa) << std::endl;

  // Get previous char for each suffix in sorted order.
  // For i == 0 wrap around and get last character.
  // for the example: "bcc$aaabb"
  auto bwt = parlay::tabulate(n+1, [&] (size_t i) -> uchar {
      auto j = sa[i];
      return (j == 0) ? ss[n] : ss[j-1];});
  
  // std::cout << parlay::map(bwt,[] (uchar x) {return (char) x;}) << std::endl;
  return bwt;
}
