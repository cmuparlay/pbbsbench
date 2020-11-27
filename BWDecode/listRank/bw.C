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
#include <set>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/internal/uninitialized_sequence.h"
#include "parlay/random.h"
#include "parlay/io.h"
#include "parlay/internal/collect_reduce.h"
#include "parlay/internal/get_time.h"
#include "bw.h"

using std::cout;
using std::endl;

// Int needs to be large enough to store s.size().
template <class Int>
ucharseq bw_decode_(ucharseq const &s) {
  parlay::internal::timer t("trans", false);
  Int n = s.size();

  struct link {Int next; uchar c;
    link(Int next, uchar c) : next(next), c(c) {}};
  
  // sort character, returning original locations in sorted order
  auto lnks = parlay::delayed_tabulate(n, [&] (size_t i) {
     return link(i, s[i]);});
  auto [links, c_] = parlay::internal::count_sort(parlay::make_slice(lnks), s, 256);
  t.next("count sort");

  // break lists into blocks
  Int block_size = 5000;

  // pick a set of about n/block_size locations as heads
  // head_flags are set to true for heads
  // links that point to a head are set to their original position + n
  // the overall first character is made to be a head
  parlay::random r(0);
  parlay::sequence<bool> head_flags(n, false);
  size_t start = links[0].next;
  head_flags[start] = true; // first char is a head
  links[0].next += n;
  parlay::parallel_for(0, n/block_size + 2, [&] (Int i) {
      size_t j = r.ith_rand(i)%n;
      auto lnk = links[j].next;
      // if not already incremented, add n (race is Ok, only inc. once)
      if (lnk < n) {
	head_flags[lnk] = true;
	links[j].next = lnk + n;
      }
    }, 1000);

  t.next("set next");

  // indices of heads;
  auto heads = parlay::pack_index<Int>(head_flags);
  t.next("pack index");

  // map over each head and follow the list until reaching the next head
  // as following the list, add characters to a buffer
  // at the end return the buffer trimmed to fit the substring exactly
  // also return pointer to the next head
  auto blocks = parlay::map(heads, [&] (Int my_head) {
      // very unlikely to take more than this much space,
      // throws an exception if it does
      Int buffer_len = block_size * 30;
      auto buffer = parlay::sequence<uchar>::uninitialized(buffer_len);
      Int i = 0;
      Int pos = my_head;
      do {
	link ln = links[pos];
	buffer[i++] = ln.c;
	if (i == buffer_len)
	  throw std::runtime_error("ran out of buffer space in bw decode");
	pos = ln.next;
      } while (pos < n);
      auto trimmed = parlay::tabulate(i, [&] (size_t j) {return buffer[j];});
      return std::make_pair(std::move(trimmed), pos % n);
    }, 1);
  t.next("follow pointers");

  // location in heads for each head in s
  auto location_in_heads = parlay::sequence<Int>::uninitialized(n);
  Int m = heads.size();
  parlay::parallel_for(0, m, [&] (Int i) {
      location_in_heads[heads[i]] = i; });
  t.next("link heads");
  
  // start at first block and follow next pointers
  // putting each into ordered_blocks
  Int pos = start;
  parlay::sequence<parlay::sequence<uchar>> ordered_blocks(m);
  for (Int i=0; i < m; i++) {
    Int j = location_in_heads[pos];
    pos = blocks[j].second;
    ordered_blocks[i] = std::move(blocks[j].first);
  }
  t.next("order heads");

  // flatten ordered blocks into final string
  auto res = parlay::flatten(ordered_blocks);
  t.next("flatten");

  // drop the first character, which is a null character
  auto rd = parlay::to_sequence(res.cut(0,res.size()-1));
  t.next("remove first");
  return rd;
}

ucharseq bw_decode(ucharseq const &s) {
  if (s.size() >= (((long) 1) << 31))
    return bw_decode_<unsigned long>(s);
  else
    return bw_decode_<unsigned int>(s);
}
