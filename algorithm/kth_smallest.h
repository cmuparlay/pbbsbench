#pragma once
#include "../parlay/primitives.h"
#include "../parlay/random.h"

namespace pbbs {
  using namespace parlay;

  template <class Seq, class Compare>
  auto kth_smallest(Seq const &s, size_t k, Compare less,
		    parlay::random r = parlay::random()) {
    using T = typename Seq::value_type;
    size_t n = s.size();
    T pivot = s[r[0]%n];
    sequence<T> smaller = filter(s, [&] (T a) {return less(a, pivot);});
    if (k < smaller.size())
      return kth_smallest(smaller, k, less, r.next());
    else {
      sequence<T> larger = filter(s, [&] (T a) {return less(pivot, a);});
      if (k >= n - larger.size())
	return kth_smallest(larger, k - n + larger.size(), less, r.next());
      else return pivot;
    }
  }

  template <class Seq, class Compare>
  auto approximate_kth_smallest(Seq const &S, size_t k, Compare less,
				parlay::random r = parlay::random()) {
    // raise exception if empty sequence?
    using T = typename Seq::value_type;
    size_t n = S.size();
    size_t num_samples = n/sqrt(n);
    auto samples = tabulate(num_samples, [&] (size_t i) -> T {
	return S[r[i]%n];});
    return parlay::sort(samples, less)[k * num_samples / n];
    //kth_smallest(samples, k * num_samples / n, less);
  }
}
