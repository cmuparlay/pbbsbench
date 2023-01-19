#pragma once
#include <parlay/parallel.h>
#include <atomic>

namespace flck {
  namespace internal {
// A dummy wrapper if no tagging is needed
template <typename V>
struct no_tagged {
  using IT = size_t;
  static inline IT init(V v) {return  (IT) v;}
  static inline V value(IT v) {return (V) v;}
  static inline IT get_tag(IT v) {return 0;}
  static bool cas(std::atomic<IT> &loc, IT oldv, V v, bool aba_free=false) {
    return loc.compare_exchange_strong(oldv, (IT) v);
  }
  static bool cas_with_same_tag(std::atomic<IT> &loc, IT oldv, V v,
				bool aba_free=false) {
    return cas(loc, oldv, v);
  }
};
  } // end namespace internal
} // end namespace flck
