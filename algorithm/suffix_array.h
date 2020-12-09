// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011-2019 Guy Blelloch, Julian Shun and the PBBS team
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

// A modified version of the Apostolico, Iliopoulos, Landau, Schieber, and Vishkin
// Suffix Array algorithm (actually originally designed as a suffix tree algorithm)
// It does O(n log n) work in the worst case, but for most inputs
//    it does O(n) work beyond a sort on fixed length keys.
// The depth is O(log^2 n) assuming the sort is within that bound.
// Each round doubles the substring lengths that have been sorted but drops any strings
// that are already sorted (i.e. have no equal strings for the current length).
// For most inputs, most strings (suffixes) drop out early.

// Supports the following interface returning a suffix array for c
//   indexT is the type of integer for the suffix indices
//   It needs to have at least log_2 n bits.  It can be unsigned.

//  // indexT is an integer type, either signed or unsigned
//  // for input of length n, but hold integers up to n+3
//  template <typename indexT>
//  parlay::sequence<indexT> suffix_array(parlay::sequence<unsigned char> const &s);

#include <math.h>
#include "../parlay/parallel.h"
#include "../parlay/primitives.h"
#include "../parlay/internal/get_time.h"

constexpr bool verbose = false;

using uchar = unsigned char;
using uint128 = unsigned __int128;

template <typename indexT>
using ipair = std::pair<indexT,indexT>;

template <typename indexT>
struct seg {
  indexT start;
  indexT length;
  seg<indexT>(indexT s, indexT l) : start(s), length(l) {}
  seg<indexT>() {}
};

template <typename indexT>
void split_segment(parlay::slice<seg<indexT>*,seg<indexT>*> segOut,
		   indexT start,
		   parlay::sequence<indexT> &ranks,
		   parlay::slice<ipair<indexT>*,ipair<indexT>*> Cs) {
  indexT l = segOut.size();
  if (l < 5000) { // sequential version

    indexT name = 0;
    ranks[Cs[0].second] = name + start + 1;
    for (indexT i=1; i < l; i++) {
      if (Cs[i-1].first != Cs[i].first) name = i;
      ranks[Cs[i].second] = name + start + 1;
    }

    name = 0;
    for (indexT i=1; i < l; i++) {
      if (Cs[i-1].first != Cs[i].first) {
	segOut[i-1] = seg<indexT>(name+start,i-name);
	name = i;
      } else segOut[i-1] = seg<indexT>(0,0);
    }
    segOut[l-1] = seg<indexT>(name+start,l-name);

  } else { // parallel version
    auto names = parlay::sequence<indexT>::uninitialized(l);

    // mark start of each segment with equal keys
    parlay::parallel_for (1, l, [&] (size_t i) {
	names[i] = (Cs[i].first != Cs[i-1].first) ? i : 0;});
    names[0] = 0;

    // scan start i across each segment ???
    parlay::scan_inclusive_inplace(names, parlay::maxm<indexT>());

    // write new rank into original location
    parlay::parallel_for (0, l, [&] (size_t i) {
	ranks[Cs[i].second] = names[i]+start+1;});

    // get starts and lengths of new segments
    parlay::parallel_for (1, l, [&] (size_t i) {
	if (names[i] == i)
	  segOut[i-1] = seg<indexT>(start+names[i-1],i-names[i-1]);
	else segOut[i-1] = seg<indexT>(0,0);
      });
    segOut[l-1] = seg<indexT>(start+names[l-1],l-names[l-1]);
  }
}

template <class indexT>
parlay::sequence<ipair<indexT>>
split_segment_top(parlay::sequence<seg<indexT>> &segOut,
		  parlay::sequence<indexT> &ranks,
		  parlay::sequence<uint128> const &Cs) {
  size_t n = segOut.size();
  auto names = parlay::sequence<indexT>::uninitialized(n);
  size_t mask = ((((size_t) 1) << 32) - 1);

  // mark start of each segment with equal keys
  parlay::parallel_for (1, n, [&] (size_t i) {
      names[i] = ((Cs[i] >> 32) != (Cs[i-1] >> 32)) ? i : 0;});
  names[0] = 0;

  // scan start i across each segment ???
  parlay::scan_inclusive_inplace(names, parlay::maxm<indexT>());

  auto C = parlay::sequence<ipair<indexT>>::uninitialized(n);
  // write new rank into original location
  parlay::parallel_for (0, n, [&] (size_t i) {
      ranks[Cs[i] & mask] = names[i]+1;
      C[i].second = Cs[i] & mask;
    });

  // get starts and lengths of new segments
  parlay::parallel_for (1, n, [&] (size_t i) {
      if (names[i] == i)
	segOut[i-1] = seg<indexT>(names[i-1],i-names[i-1]);
      else segOut[i-1] = seg<indexT>(0,0);
    });
  segOut[n-1] = seg<indexT>(names[n-1],n-names[n-1]);

  return C;
}

template <class indexT, class UCharRange>
parlay::sequence<indexT> suffix_array(UCharRange const &ss) {
  parlay::internal::timer sa_timer("Suffix Array", false);
  size_t n = ss.size();

  // renumber characters densely
  // start numbering at 1 leaving 0 to indicate end-of-string
  size_t pad = 48;
  parlay::sequence<indexT> flags(256, (indexT) 0);
  parlay::parallel_for (0, n, [&] (size_t i) {
      if (!flags[ss[i]]) flags[ss[i]] = 1;}, 1000);
  auto add = [&] (indexT a, indexT b) {return a + b;};
  indexT m;
  std::tie(flags, m) = parlay::scan(flags, parlay::make_monoid(add,(indexT) 1));

  // pad the end of string with 0s
  auto s = parlay::tabulate(n + pad, [&] (size_t i) -> uchar {
      return (i < n) ? flags[ss[i]] : 0;});

  if (verbose) std::cout << "distinct characters = " << m-1 << std::endl;

  // pack characters into 128-bit word, along with the location i
  // 96 bits for characters, and 32 for location
  double logm = log2((double) m);
  indexT nchars = floor(96.0/logm);

  auto Cl = parlay::tabulate(n, [&] (size_t i) -> uint128 {
      uint128 r = s[i];
      for (indexT j=1; j < nchars; j++) r = r*m + s[i+j];
      return (r << 32) + i;
    });
  sa_timer.next("copy into 128bit int");

  // sort based on packed words ??
  parlay::sort_inplace(Cl, std::less<uint128>());
  sa_timer.next("sort");

  // identify segments of equal values
  auto ranks = parlay::sequence<indexT>::uninitialized(n);
  auto seg_outs = parlay::sequence<seg<indexT>>::uninitialized(n); 
  parlay::sequence<ipair<indexT>> C = split_segment_top(seg_outs, ranks, Cl);
  Cl.clear();
  sa_timer.next("split top");

  indexT offset = nchars;
  uint round =0;
  indexT nKeys = n;

  // offset is how many characters for each suffix have already been sorted
  // each round doubles offset so there should be at most log n rounds
  // The segments keep regions that have not yet been fully sorted
  while (1) {
    if (round++ > 40)
      throw std::runtime_error("Suffix Array: internal error, too many rounds");

    auto is_seg = [&] (seg<indexT> s) {return s.length > 1;};
    // only keep segments that are longer than 1 (otherwise already sorted)
    parlay::sequence<seg<indexT>> Segs = parlay::filter(seg_outs.cut(0,nKeys), is_seg);
    indexT nSegs = Segs.size();
    if (nSegs == 0) break;
    sa_timer.next("filter and scan");

    auto offsets = parlay::sequence<indexT>::uninitialized(nSegs);
    parlay::parallel_for (0, nSegs, [&] (size_t i) {
	indexT start = Segs[i].start;
	indexT l = Segs[i].length;
	auto Ci = C.cut(start, start + l);
	offsets[i] = l;

	// grab rank from offset locations ahead
	parlay::parallel_for (0, l, [&] (size_t j) {
	    indexT o = Ci[j].second + offset;
	    Ci[j].first = (o >= n) ? 0 : ranks[o];
	  }, 100);

	// sort within each segment based on ranks
	auto less = [&] (ipair<indexT> A, ipair<indexT> B) {
	  return A.first < B.first;};
	parlay::sort_inplace(Ci, less);
      });
    sa_timer.next("sort");

    // starting offset for each segment ???
    nKeys = parlay::scan_inplace(offsets, parlay::addm<indexT>());

    // Split each segment into subsegments if neighbors differ.
    parlay::parallel_for (0, nSegs, [&] (size_t i) {
	indexT start = Segs[i].start;
	indexT l = Segs[i].length;
	indexT o = offsets[i];
	split_segment(seg_outs.cut(o, o + l),
		      start,
		      ranks,
		      C.cut(start, start+l));
      }, 100);
    sa_timer.next("split");
    
    if (verbose)
      std::cout << "length: " << offset << " keys remaining: " << nKeys << std::endl;
    
    offset = 2 * offset;
  }
  parlay::parallel_for (0, n, [&] (size_t i) {
      ranks[i] = C[i].second;});
  return ranks;
}
