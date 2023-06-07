#pragma once

// The flock library.

// Public interface for lock:
//   with_lock(thunk_returning_val) -> val
//   try_lock(thunk_returning_boolean) -> boolean
//   try_lock_result(thunk_returning_val) -> optional<val>
//   wait_lock() -> void
//   is_locked() -> bool

// Supported data types

//  atomic<T> :
//    atomic(v) : constructor
//    atomic() : default constructor
//    load() : return value
//    store(v) : store value
//    cam(old_v, new_v) : cas but does not return whether successful
//    ** the following two are the same as load and store if used
//       with regular locks
//    read() : return value, but not idempotent
//    init(v) : like store(v) but only used before anyone other thread
//              has a handle to this mutable
//    ** the following only needed for snapshotting and only
//       if trying to optimized code.
//    validate() : ensures it has a timestamp

//  ** The following is if a value is only changed once from its
//     initial value
//  write_once<T> :
//    write_once(v) : constructor
//    write_once() : default constructor
//    load() : return value
//    store(v) : store value
//    init(v) : like store(v) but only used before anyone other thread
//              has a handle to this mutable

//  memory_pool<T> :
//    memory_pool() : constructor
//    new_obj(constructor args) -> *T : new object of type T
//    retire(*T) -> void : returns memory to pool
//    ** the following only needed for lock_free locks
//    new_init(f, constructor args) : applies f to constructed object
//    ** Statistics and others (can be noops)
//    stats()
//    shuffle()
//    reserve()
//    clear()


#ifdef NoHelp  // use regular locks
#include "spin_lock.h"
#include "lock_types.h"
#else  // use lock free locks
#include "lf_lock.h"
#include "lf_types.h"
#endif

namespace flck {

// This selects between using a hashlock or inline lock at compile time
// For the hashlock the address of the structure is hashed to one of a
// fixed number of lock locations and there is no need for a lock per
// structure.  However, the hashing can create lock cycles so this
// cannot be used with a strict lock, just a try_lock.

#ifdef HashLock
struct lock {
private:
  static const int bucket_bits = 16;
  static const size_t mask = ((1ul) << bucket_bits) - 1;
  static std::vector<internal::lock> locks;
  internal::lock* get_lock() {
    return &locks[parlay::hash64_2((size_t) this) & mask];}
public:
  template <typename F>
  bool try_lock(F f) {
    return get_lock()->try_lock(f);}
  void wait_lock() { get_lock()->wait_lock(); }
  bool is_locked() { return get_lock()->is_locked(); }
};

std::vector<internal::lock> lock::locks{1ul << bucket_bits};

#else 
 using lock = internal::lock;
#endif
  thread_local long try_time_taken = -1;

  template <typename F>
  auto try_loop(const F& f, int delay = 200, const int max_multiplier = 10) {
    int multiplier = 1;
    int cnt = 0;
    while (true)  {
      if (cnt++ == 1000000000/(delay*max_multiplier)) {
	std::cout << "problably in an infinite retry loop" << std::endl;
	abort(); 
	//deadlock = true;
      }
      auto r = f();
      //if (try_time_taken != -1) delay = try_time_taken;
      if (r.has_value()) return *r;
      multiplier = std::min(2*multiplier, max_multiplier);
      for (volatile int i=0; i < delay * multiplier; i++);
    }
  }

  template <typename F>
  auto try_loop_(int rounds, const F& f, int delay = 200, const int max_multiplier = 10) {
    int multiplier = 1;
    int cnt = 0;
    while (cnt++ < rounds)  {
      auto r = f();
      if (r.has_value()) return r;
      multiplier = std::min(2*multiplier, max_multiplier);
      for (volatile int i=0; i < delay * multiplier; i++);
    }
    using RT = decltype(f());
    return RT();
  }

} // namespace flck
