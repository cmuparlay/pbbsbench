#pragma once
// Public interface for lock:
//   with_lock(thunk_returning_val) -> val
//   try_lock(thunk_returning_boolean) -> boolean
//   try_lock_result(thunk_returning_val) -> optional<val>
//   wait_lock() -> void
//   is_locked() -> bool
//   is_self_locked() -> bool

#include <atomic>
#include <functional>
#include <assert.h>
#include "lf_log.h"
#include "lf_types.h"
#include "acquired_pool.h"

namespace flck {
namespace internal {

// used for reentrant locks
  static thread_local size_t current_id = 1000000;
  inline size_t get_current_id() {
    if (current_id == 1000000)
      current_id = epoch::internal::worker_id();
    return current_id;
  }
  inline void set_current_id(size_t id) {
    current_id = id;
  }

// user facing lock
using lock_entry_ = size_t;
struct lock; // for back reference

// stores the thunk along with the log
struct descriptor {
  std::function<void()> f;  // the thunk to run
  bool done;  // set when done
  bool freed; // just for debugging
  // Used for memory management to indicate the thunk is being helped.
  // It lives beyond the lifetime of the descriptor.
  // In particular a helper might set this even if being used reused.
  // This is OK since setting to true is always semantically OK.
  // Being false just allows more efficient (i.e. immediate) reclamation.
  std::atomic<bool> acquired;   
  int thread_id; // used to detect reentrant locks, inherited by helpers
  lock_entry_ current; // currently not used
  std::atomic<size_t> counter; // used for aba-prevention when taking a lock
  long epoch_num; // the epoch when initially created, inherited by helpers
  log_array lg_array; // the log itself
  
  template <typename Thunk>
  descriptor(Thunk& g, long counter) : 
    f([=] {g();}), done(false), freed(false), counter(counter) {
    lg_array.init();
    epoch_num = epoch::internal::get_epoch().get_my_epoch();
    thread_id = get_current_id();
  }

  ~descriptor() { // just for debugging
    assert(!freed);
    freed = true;
  }

  void operator () () {
    assert(!freed);
    // run f using log based on lg_array
    with_log(Log(&lg_array,0), f);
    done = true;
    //std::atomic_thread_fence(std::memory_order_seq_cst);
  }
};

#ifdef EpochDescriptor
// use epoch based collector to reclaim descriptors
memory_pool<descriptor> descriptor_pool;

#else
// Use optimized collector to reclaim descriptors.
// If noone is helping (acquire flag is false), then reclaim
// immediately, otherwise send to epoch-based collector.
// When helping the acquired flag needs to be set.

inline memory_pool<descriptor,acquired_pool<descriptor>>& get_descriptor_pool() {
  static memory_pool<descriptor,acquired_pool<descriptor>> descriptor_pool;
  return descriptor_pool;
}

#endif

struct lock {
public:
    using lock_entry = lock_entry_;
private:
  // each lock entry will be a pointer to a descriptor when locked
  // and a counter when unlocked. The counter is to prevent ABA problems.
  std::atomic<lock_entry> lck;
  lock_entry load() {return lg.commit_value(lck.load()).first;}
  lock_entry read() {return lck.load();}

  // used to take lock for version with helping
  bool cas(lock_entry oldl, descriptor* d) {
    lock_entry current = read();
    if (current != oldl) return false;
    return lck.compare_exchange_strong(oldl, (lock_entry) d);
  }
  
  void clear(descriptor* d) {
    lock_entry current = lck.load();
    if (current == (lock_entry) d)
      // true indicates this is ABA free since current cannot
      // be reused until all helpers are done with it.
      lck.compare_exchange_strong(current, d->counter+1);
  }

  bool is_locked_(lock_entry le) { return (le & (1ull<<63)) == 0;}
  descriptor* remove_tag(lock_entry le) { return (descriptor*) le;}
  bool lock_is_self(lock_entry le) {
    return ((int)get_current_id()) == remove_tag(le)->thread_id;
  }

  // runs thunk in appropriate epoch and after it is acquired
  bool help_descriptor(lock_entry le, bool recursive_help=false) {
    if (!recursive_help && helping) return false;
    descriptor* desc = remove_tag(le);
    bool still_locked = (read() == le);
    if (!still_locked) return false;
    long my_epoch = epoch::internal::get_epoch().get_my_epoch();
    long other_epoch = desc->epoch_num;
    if (other_epoch < my_epoch)
      epoch::internal::get_epoch().set_my_epoch(other_epoch); // inherit epoch of helpee
    int my_id = get_current_id(); 
    set_current_id(desc->thread_id);   // inherit thread id of helpee
    get_descriptor_pool().acquire(desc);  // mark descriptor as acquired
    still_locked = (read() == le);
    if (still_locked) {
      bool hold_h = helping; 
      helping = true; // mark as in helping mode
      (*desc)();      // run thunk to be helped
      clear(desc);    // unset the lock
      helping = hold_h; // reset helping mode
    }
    set_current_id(my_id); // reset thread id
    epoch::internal::get_epoch().set_my_epoch(my_epoch); // reset to my epoch
    return still_locked; // return true if did helping
  }

public:
  lock() : lck((1ull<<63)) {} // initialize in unlocked state
  
  // waits until current owner of lock releases it (unless self owned)
  void wait_lock() {
    lock_entry current = load();
    if (!is_locked_(current) || lock_is_self(current)) return;
    // last arg needs to be true, otherwise might not help and clear
    help_descriptor(current, true);
  }

  bool is_locked() { return is_locked_(load());}
  bool is_self_locked() {
    lock_entry current = load();
    return is_locked_(current) && lock_is_self(current);
  }
  lock_entry lock_load() {return load();}
  bool unchanged(lock_entry le) {return le == load();}

  
  // This is safe to be used inside of another lock (i.e. it is
  // idempotent, kind of).  The key components to making it effectively
  // idempotent is using an idempotent new_obj, and checking if it is
  // done (my_thunk->done).  It is lock free if no cycles in lock
  // ordering, and otherwise can deadlock.
  template <typename Thunk>
  auto with_lock(Thunk f) {
    using RT = decltype(f());
    static_assert(sizeof(RT) <= 4 || std::is_pointer<RT>::value,
		  "Result of with_lock must be a pointer or at most 4 bytes");
    lock_entry current = read();
    bool locked = is_locked_(current);

    // idempotently allocate descriptor
    auto [my_descriptor, i_own] = get_descriptor_pool().new_obj_acquired(f, locked ? 0 : current);
    
    // if already retired, then done
    if (get_descriptor_pool().is_done(my_descriptor)) {
      auto ret_val = get_descriptor_pool().done_val_result<RT>(my_descriptor);
      assert(ret_val.has_value()); // with_lock is guaranteed to succeed
      return ret_val.value(); 
    }
    
    
    while (true) {
      size_t old_count = my_descriptor->counter;
      if(!locked && old_count < current) {
        if(!my_descriptor->counter.compare_exchange_strong(old_count, current)) {
          current = old_count;
          locked = is_locked_(current);
        }
      }
      if (my_descriptor->done // already done
          || remove_tag(current) == my_descriptor // already acquired
          || (!locked && cas(current, my_descriptor))) { // try to acquire

        // run the body f with the log from my_descriptor
        RT result = with_log(Log(&my_descriptor->lg_array,0), [&] {return f();});

        // mark as done and clear the lock
        my_descriptor->done = true;
        clear(my_descriptor);

        // retire the descriptor saving the result in the enclosing
        // descriptor, if any
        get_descriptor_pool().retire_acquired_result(my_descriptor, i_own,
                                               std::optional<RT>(result));
        return result;
      } else if (locked) {
        help_descriptor(current);
      }
      current = read();
      locked = is_locked_(current);
    }
  }

  void unlock() {
    lock_entry current = load();
    lck.compare_exchange_strong(current, ((descriptor*) current)->counter+1);
  }

  // The thunk returns a value
  // The try_lock_result returns an optional value, which is empty if it fails
  template <typename Thunk>
  auto try_lock_result(Thunk f, bool do_help=true) {
    using RT = decltype(f());
    std::optional<RT> result = {};
    lock_entry current = load();

    // check if reentrant lock (already locked by self)
    if (is_locked_(current) && lock_is_self(current)) 
      return std::optional(f()); // if so, run without acquiring

    // Idempotent allocation of descriptor.  
    auto [my_descriptor, i_own] = get_descriptor_pool().new_obj_acquired(f, is_locked_(current) ? 0 : current);
    
    // if descriptor is already retired, then done and return value
    if (get_descriptor_pool().is_done(my_descriptor)) 
      return get_descriptor_pool().done_val_result<RT>(my_descriptor);
        
    if (!is_locked_(current)) {
      // use a CAS to try to acquire the lock
      cas(current, my_descriptor);

      // This could be a load() without the my_descriptor->done test.
      // Using read() is an optimization to avoid a logging event.
      current = read();
      if (my_descriptor->done || remove_tag(current) == my_descriptor) {

        // run f with log from my_descriptor
        result = with_log(Log(&my_descriptor->lg_array,0),
                          [&] {return f();});

        // mark as done and clear the lock
        my_descriptor->done = true;
        clear(my_descriptor);
      }
    } else {
      if (do_help) help_descriptor(current);
    }

    // retire the thunk
    get_descriptor_pool().retire_acquired_result(my_descriptor, i_own, result);
    return result;
  }

  // A wrapper that returns true if the try_lock succeeds and its thunk
  // returns true.
  template <typename Thunk>
  auto try_lock(Thunk f, bool do_help=true) {
    auto result = try_lock_result(f, do_help);
    return result.has_value() && result.value();
  }

  // for compatibility with transactions, otherwise do not use
  template <typename Thunk>
  auto read_lock(Thunk f) {return f();}

};

} // namespace internal
} // namespace flck
