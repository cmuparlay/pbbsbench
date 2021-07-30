#include "parlay/primitives.h"
#include <unordered_set>

template <class T>
parlay::sequence<T> dedup(parlay::sequence<T> const &A) {
  std::unordered_set<T> deduped(A.size());
  for (size_t i=0; i<A.size(); i++) {
    deduped.insert(A[i]);
  }
  parlay::sequence<T> result;
  result.reserve(deduped.size());
  for (const auto& val : deduped) {
    result.push_back(val);
  }
  return result;
}
