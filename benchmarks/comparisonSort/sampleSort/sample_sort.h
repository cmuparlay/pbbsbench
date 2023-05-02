#include <algorithm>
#include <functional>
#include <random>

#include <parlay/parallel.h>
#include <parlay/primitives.h>
#include <parlay/random.h>
#include <parlay/internal/get_time.h>
#include <parlay/slice.h>
#include <parlay/utilities.h>
#include <parlay/internal/uninitialized_sequence.h>

#include "heap_tree.h"

// **************************************************************
// Sample sort
// A generalization of quicksort to many pivots.
// This code picks up to 256 pivots by randomly selecting and
// then sorting them.
// It then puts the keys into buckets depending on which pivots
// they fall between and recursively sorts within the buckets.
// Makes use of a parlaylib bucket sort for the bucketing.
//
// This is taken from parlaylib/examples/samplesort.h with some optimizations.
// In particular this uses uninitialized sequences and relocation, which
// avoids copying the keys.  Useful for strings.
// Also uses parlay's internal quicksort (faster for strings)
// **************************************************************

template <typename assignment_tag, typename Range, typename Less>
void sample_sort_(Range in, Range out, Less less, bool stable=false, int level=1) {
  long n = in.size();
  parlay::internal::timer t("sample", level==1);
  using T = typename Range::value_type;
  using bucket_key_t = unsigned short;
  
  // for the base case (small or recursion level greater than 2)
  long cutoff = 1024;
  if (n <= cutoff || (level > 2)) { // && n <= 1 << 17)) {
    //return;
    if (in.begin() != out.begin())
      parlay::uninitialized_relocate_n(out.begin(), in.begin(), n);
    if (stable)
      std::stable_sort(out.begin(), out.end(), less);
    else
      if (sizeof(T) > 8)
	parlay::internal::quicksort(out.begin(), n, less);
      else std::sort(out.begin(), out.end(), less);
    return;
  }

  // number of bits in bucket count (e.g. 8 would mean 256 buckets)
  int max_bits = (level == 1) ? 10 : 8;
  int bits = std::min<long>(max_bits, parlay::log2_up(n)-parlay::log2_up(cutoff)+1);
  long num_buckets = 1 << bits;

  // over-sampling ratio: keeps the buckets more balanced
  int over_ratio = 8;

  // create an over sample and sort it using std:sort
  parlay::random_generator gen;
  std::uniform_int_distribution<long> dis(0, n-1);
  auto oversample = parlay::tabulate(num_buckets * over_ratio, [&] (long i) {
    auto r = gen[i];
    return in[dis(r)];}, 1000);
  std::sort(oversample.begin(), oversample.end(), less);

  // sub sample to pick final pivots (num_buckets - 1 of them)
  auto pivots = parlay::tabulate(num_buckets-1, [&] (long i) {
						  return oversample[(i+1)*over_ratio];}, 1000);

  // check if any duplicates among the pivots
  bool duplicates = false;
  for (int i=0; i < num_buckets-2; i++)
    if (!less(pivots[i],pivots[i+1])) duplicates = true;
  
  // put pivots into efficient search tree and find buckets id for the input keys
  heap_tree ss(pivots);
  parlay::sequence<bucket_key_t> bucket_ids;
  if (duplicates)
    // if duplicates put keys equal to a pivot in next bucket
    // this ensures all keys equaling a duplicate are in a bucket by themselves
    bucket_ids = parlay::tabulate(n, [&] (long i) -> bucket_key_t {
			  auto b = ss.find(in[i], less);
			  return b + ((b < num_buckets-1) && !less(in[i],pivots[b]));}, 1000);
  else
    bucket_ids = parlay::tabulate(n, [&] (long i) -> bucket_key_t {
				       return ss.find(in[i], less);}, 1000);
  t.next("bucket id");
   
  // sort into the buckets
  auto keys = parlay::internal::uninitialized_sequence<T>(n);
  auto keys_slice = parlay::make_slice(keys);
  auto offsets = parlay::internal::count_sort<assignment_tag>(in, keys_slice,
							      parlay::make_slice(bucket_ids), num_buckets).first;
  t.next("count sort");
  
  // now recursively sort each bucket
  parlay::parallel_for(0, num_buckets, [&, &keys = keys, &offsets = offsets] (long i) {
    long first = offsets[i];
    long last = offsets[i+1];
    if (last-first == 0) return; // empty
    // if duplicate keys among pivots don't sort all-equal buckets
    if (i == 0 || i == num_buckets - 1 || less(pivots[i-1], pivots[i]))
      sample_sort_<parlay::uninitialized_relocate_tag>(keys_slice.cut(first,last), out.cut(first,last),
						       less, stable, level+1);
    else parlay::uninitialized_relocate_n(out.begin()+first, keys_slice.begin()+first, last-first);
  }, 1);
}

// A version that returns sequence in the input.
// It does use O(n) temporary memory.
template <typename Range, typename Less = std::less<>>
void sample_sort_inplace(Range& in, Less less = {}) {
  auto ins = parlay::make_slice(in);
  sample_sort_<parlay::uninitialized_relocate_tag>(ins, ins, less);
}

//  A version that does not destroy the input
template <typename Range, typename Less = std::less<>>
auto sample_sort(const Range& in, Less less = {}) {
  using T = typename Range::value_type;
  auto ins = parlay::make_slice(in);
  auto outs = parlay::sequence<T>::uninitialized(ins.size());
  sample_sort_<parlay::uninitialized_copy_tag>(ins, parlay::make_slice(outs), less);
  return outs;
}

//  A version that does not destroy the input
template <typename Range, typename Less = std::less<>>
auto stable_sort(const Range& in, Less less = {}) {
  using T = typename Range::value_type;
  auto ins = parlay::make_slice(in);
  auto outs = parlay::sequence<T>::uninitialized(ins.size());
  sample_sort_<parlay::uninitialized_copy_tag>(ins, parlay::make_slice(outs), less, true);
  return outs;
}
