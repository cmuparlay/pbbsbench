// This selects between using versioned objects or regular objects
// Versioned objects are implemented as described in:
//   Wei, Ben-David, Blelloch, Fatourou, Rupert and Sun,
//   Constant-Time Snapshots with Applications to Concurrent Data Structures
//   PPoPP 2021
// They support snapshotting via version chains, and without
// indirection, but pointers to objects (ptr_type) must be "recorded
// once" as described in the paper.
#pragma once
#include <parlay/parallel.h>
#include <parlay/sequence.h>
#include "../flock/flock.h"

namespace verlib {
 bool strict_lock = false;
}

#ifdef Versioned

// versioned objects, ptr_type includes version chains
#ifdef Recorded_Once
#include "versioned_recorded_once.h"
#elif FullyIndirect
#include "versioned_indirect.h"
#elif Transactional
#include "../transactions/transactions.h"
#elif GenSnapshot
#include "versioned_generalized.h"
#else
#include "versioned_hybrid.h"
#endif

# else // Not Versioned

namespace verlib {

  struct versioned {};

  template <typename T>
  using versioned_ptr = flck::atomic<T*>;
  using atomic_bool = flck::atomic<bool>;
  using flck::lock;
  using flck::memory_pool;
  template <typename A, typename B, typename C>
  bool validate(const A& a, const B& b, const C& c) {return true;}

  template <typename F>
  auto with_snapshot(F f, bool unused_parameter=false) {
    return flck::with_epoch([&] { return f();});
  }
  template <typename F>
  auto do_now(F f) { return f();}  
}

#endif

namespace verlib {
  using flck::with_epoch;
}
