#pragma once
#include <atomic>
#include <vector>
#include <parlay/alloc.h>
#include "epoch.h"

// A trivial wrapper around a memory pool so that retire will call
// "destruct" instead of "retire" if the object's "acquired" field is
// false.  The object type xT must have an atomic boolean field called
// acquired.  Note that "destruct" frees the memory immediately while
// "retire" puts it aside and does not free until it is safe to do so
// (e.g. two epochs later in an epoch-based pool).  Object needs
// acquired field and if acquired after it is retired user then needs
// to check by some means that not retired.  Used to manage
// descriptors in the lock-free locks.  Any helpers will acquire the
// descriptor.  They first acquire (possibly after the descriptor is
// retired) and then check that the lock still points to the
// descriptor.  On retiring a descriptor it is first removed from the
// lock, then retired.  Note this is conservative in that acquired can
// be accidentally true (on initialization, or due to the race
// described above).
namespace flck {
  namespace internal {
template <typename xT>
struct acquired_pool {
  using T = xT;
  mem_pool<T> pool;
  void reserve(size_t n) { pool.reserve(n);}
  void shuffle(size_t n) { pool.shuffle(n);}
  void stats() { pool.stats();}
  void clear() { pool.clear();}
  void destruct(T* ptr) { pool.destruct(ptr); }
  void acquire(T* ptr) {
    if (ptr->acquired == false)
      ptr->acquired = true; }
  template <typename ... Args>
  T* new_obj(Args... args) { return pool.new_obj(args...);} 

  void retire(T* ptr) {
    bool x = (ptr->acquired).load();
    if (!x) destruct(ptr);
    else {
      ptr->acquired = false;
      pool.retire(ptr);
    }
  }
};
  } // namespace internal
} // namespace flck
