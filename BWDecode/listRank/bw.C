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
#include "parlay/internal/uninitialized_sequence.h"
#include "parlay/random.h"
#include "parlay/io.h"
#include "parlay/internal/collect_reduce.h"
#include "common/get_time.h"
#include "bw.h"

using std::cout;
using std::endl;

template <class Int>
ucharseq bw_decode_(ucharseq const &s, size_t head) {
  timer t("trans", true);

  Int n = s.size();
  ucharseq ss = s;
  t.next("copy");

  // Set head location to 0 temporarily so when sorted it ends up
  // first this is the right position since the last character in the
  // original string is the first when sorted by suffix, and head is
  // its "suffix" with wraparound
  uchar a = ss[head];
  ss[head] = 0;

  // sort character, returning original locations in sorted order
  auto [links, counts] = parlay::internal::count_sort(parlay::make_slice(parlay::iota<Int>(n)), ss, 256);
  ss[head] = a;
  t.next("count sort");

  Int block_bits = 10;
  Int block_size = 1 << block_bits;
  parlay::random r(0);
  parlay::sequence<bool> head_flags(n, false);
  t.next("head flags");
  
  head_flags[head] = true;  // the overall head is made to be a head
  links[0] += n;   // this points to the original head

  // pick a set of n/block_size locations as heads
  // links that point to a head are set to their original position + n
  parlay::parallel_for(0, n/block_size + 2, [&] (Int i) {
      size_t j = r.ith_rand(i)%n;
      if (!head_flags[links[j]]) {
	head_flags[links[j]] = true;
	links[j] += n;
      }});

  t.next("set next");

  // indices of heads;
  auto heads = parlay::pack_index<Int>(head_flags);
  t.next("pack index");

  // map over each head and follow the list until reaching the next head
  // as following the list, add characters to a buffer
  // at the end return the buffer trimmed to fit the substring exactly
  // also return pointer to the next head
  auto blocks = parlay::map(heads, [&] (Int my_head) {
      // very unlikely to take more than this much space, throws an exception if it does
      Int buffer_len = block_size * 30;
      auto buffer = parlay::sequence<uchar>::uninitialized(buffer_len);
      Int i = 0;
      Int pos = my_head;
      do {
	buffer[i++] = ss[pos];
	if (i == buffer_len)
	  throw std::runtime_error("ran out of buffer space in bw decode");
	pos = links[pos];
      } while (pos < n);
      auto trimmed = parlay::tabulate(i, [&] (size_t j) {return buffer[j];});
      return std::make_pair(std::move(trimmed), pos - n);
    });
  t.next("follow pointers");

  // location in heads for each head in s
  auto location_in_heads = parlay::sequence<Int>::uninitialized(n);
  Int m = heads.size();
  parlay::parallel_for(0, m, [&] (Int i) {
      location_in_heads[heads[i]] = i; });
  t.next("link heads");

  // start at first block and follow next pointers
  // putting each into ordered_blocks
  Int pos = head;
  parlay::sequence<parlay::sequence<uchar>> ordered_blocks(m);
  for (Int i=0; i < m; i++) {
    Int j = location_in_heads[pos];
    ordered_blocks[i] = std::move(blocks[j].first);
    pos = blocks[j].second;
  }
  t.next("order heads");

  // flatten ordered blocks into final string
  auto res = parlay::flatten(ordered_blocks);
  t.next("flatten");
  return res;
}

ucharseq bw_decode(ucharseq const &s, size_t head) {
  if (s.size() >= (((long) 1) << 31))
    return bw_decode_<unsigned long>(s, head);
  else
    return bw_decode_<unsigned int>(s, head);
}
