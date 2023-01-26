#pragma once
// Public interface for lock:
//   with_lock(thunk_returning_val) -> val
//   try_lock(thunk_returning_boolean) -> boolean
//   try_lock_result(thunk_returning_val) -> optional<val>
//   wait_lock() -> void
//   is_locked() -> bool

#include <parlay/parallel.h> // needed for worker_id

#include <atomic>
#include<chrono>
#include<thread>
#include<optional>
#include<sstream>

namespace flck {
  namespace internal {
    
// used for reentrant locks
static thread_local size_t current_id = parlay::worker_id();

// This lock keeps track of how many times it was taken and who has it
struct lock {
public:

  // The low 32 bits are a counter of how many times lock has been
  // taken. This is mostly meant for debugging and currently not used
  // other than to indicate whether locked or not.  An odd number
  // means it is locked, and an even unlocked.  Bits [32-48) are used
  // to store one more than the thread id of who has the lock.  This
  // is to identify and allow self locking.
  struct lock_entry {
    size_t le;
    lock_entry(size_t le) : le(le) {}
    lock_entry() : le(0) {}
    bool is_locked() { return (le % 2 == 1);}
    size_t get_count() {return le & ((1ul << 32) - 1);}
    lock_entry take_lock() {
      return lock_entry(((current_id + 1ul) << 32) | get_count()+1);}
    lock_entry release_lock() { return lock_entry(get_count()+1);}
    size_t get_procid() { return (le >> 32) & ((1ul << 16) - 1);}
    bool is_self_locked() { return current_id + 1 == get_procid();}
  };

private:
  std::atomic<lock_entry> lck;
public:
  lock() : lck(lock_entry()) {}

  bool is_locked() { return lck.load().is_locked();}
  bool is_self_locked() { return lck.load().is_self_locked();}
  lock_entry lock_load() {return lck.load();}
  bool unchanged(lock_entry le) {return le.le == lck.load().le;}
  
  void wait_lock() {
    lock_entry current = lck.load();
    while (current.is_locked() && !current.is_self_locked())
      current = lck.load();
  }

  template <typename Thunk>
  auto try_lock_result(Thunk f) {
    using RT = decltype(f());
    lock_entry current = lck.load();
    if (!current.is_locked()) { // unlocked
      lock_entry newl = current.take_lock();
      if (lck.compare_exchange_strong(current, newl)) {
	RT result = f();
	lck = newl.release_lock();  // release lock
	return std::optional<RT>(result); 
      } else {
	return std::optional<RT>(); // fail
      }
    } else if (current.is_self_locked()) {// reentry
      return std::optional<RT>(f());
    } else {
      return std::optional<RT>(); // fail
    }
  }

  template <typename Thunk>
  auto try_lock(Thunk f) {
    auto result = try_lock_result(f);
    return result.has_value() && result.value();
  }
    
  template <typename Thunk>
  auto with_lock(Thunk f) {
    int retrys = 0;
    while(true) {
      auto result = try_lock_result(f);
      if(result.has_value()) return result.value();
      retrys++;
      // if(retrys % 100 == 0) {
        // std::stringstream ss;
        // ss << "thread " << current_id << " locking: " << &lck << " attempt: " << retrys << std::endl;
        // std::cerr << ss.str();
      // }
    }
  }
};
  } // namespace  internal
} //namespace flck
