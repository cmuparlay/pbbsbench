#pragma once
#include <limits>
#include <atomic>
#include "../flock/epoch.h"
#include <x86intrin.h>

namespace verlib {
  using TS = long;

thread_local int read_delay = 1;
thread_local int write_delay = 1;

  TS rdtsc(){
      unsigned int lo,hi;
      __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
      return (TS) (((uint64_t)(hi & ~(1<<31)) << 32) | lo);
  }

  thread_local TS last_ts = 0;

  // Use TL2 stamp
struct alignas(128) timestamp_tl2 {
  std::atomic<TS> stamp;
  int delay;

  // bottom thread_bits bits are used for the thread id;
  static constexpr int thread_bits = 10;
  static constexpr TS step = (1l << thread_bits);
  static constexpr TS mask = step - 1;

  TS get_stamp() {return stamp.load();}

  // mask out the processor id then compare
  bool less(TS a, TS b) {
    return ((a >> thread_bits) < (b >> thread_bits));
  }
  
  // check on lazy read if stamp is ok -- if "equal" then abort
  TS equal(TS a, TS b) {
    return ((a >> thread_bits) == (b >> thread_bits) && a != b);
  }

  // should not get here if speculative
  TS get_read_stamp() {
    std::abort(); 
    TS ts = stamp.load();
    increment_stamp(ts);
    return ts;
  }

  // apply if read_only is aborted or as part of get_write_stamp
  TS increment_stamp(TS ts) {
    TS id = epoch::internal::worker_id();
    TS new_ts = ((ts + step) & ~mask) | id;
    for (volatile int i = 1; i < delay; i++) {}
    if (ts == stamp.load()) 
      stamp.compare_exchange_strong(ts, new_ts);
    last_ts = new_ts;
    return new_ts;
  }

  // applied on every write
  TS get_write_stamp() {
    TS ts = stamp.load();
    if (ts == last_ts) return increment_stamp(ts);
    TS id = epoch::internal::worker_id();
    TS new_ts = (ts & ~mask) | id;
    last_ts = new_ts;
    return new_ts;
  }

  timestamp_tl2(int d = 20) : stamp(step), delay(d) {}
};

  // Use hardware clock with update ensured on read
struct alignas(64) timestamp_read_hw {
  TS get_stamp() {return rdtsc();}

  bool less(TS a, TS b) { return a < b;}
  bool equal(TS a, TS b) { return a == b;}
  
  
  TS get_read_stamp() {
    TS ts = rdtsc();
    while(ts == rdtsc()) {}
    std::atomic_thread_fence(std::memory_order_seq_cst);
    return ts;
  }
  TS get_write_stamp() {
    std::atomic_thread_fence(std::memory_order_seq_cst);
    return rdtsc();}
  timestamp_read_hw() {}
};

  // Use hardware clock with update ensured on write
struct alignas(64) timestamp_write_hw {
  TS get_stamp() {return rdtsc();}
  bool less(TS a, TS b) { return a < b;}
  bool equal(TS a, TS b) { return a == b;}
  TS get_read_stamp() {
    std::atomic_thread_fence(std::memory_order_seq_cst);
    return rdtsc();}
  TS get_write_stamp() {
    TS ts = rdtsc();
    while(ts == rdtsc()) {}
    std::atomic_thread_fence(std::memory_order_seq_cst);
    return ts+1;
  }
  timestamp_write_hw() {}
};


struct alignas(128) timestamp_read {
  std::atomic<TS> stamp;
  alignas(128) int delay;
  bool less(TS a, TS b) { return a < b;}
  bool equal(TS a, TS b) { return a == b;}
  
  TS get_stamp() {return stamp.load();}

  TS get_read_stamp() {
    TS ts = stamp.load();
    increment_stamp(ts);
    return ts;
  }

  void increment_stamp(TS ts) {
    if(delay == -1) {
      // delay to reduce contention
      for(volatile int i = 1; i < read_delay; i++) {}
      // std::atomic_thread_fence(std::memory_order_seq_cst);

      // only update timestamp if has not changed
      if (stamp.load() == ts) {
        if(stamp.fetch_add(1) == ts) {
	  //if (stamp.compare_exchange_strong(ts,ts+1)) {
          if(read_delay >= 2) read_delay /= 2;
        }
        else {
          if(read_delay < 256) read_delay *= 2;
        }
      }
    } else {
      for(volatile int i = 1; i < delay; i++) {}
      if (stamp.load() == ts)
        stamp.fetch_add(1);
    }
  }

  TS get_write_stamp() {return stamp.load();}
  timestamp_read(int d = -1) : stamp(1) {
    auto cstr = std::getenv("READ_DELAY");
    if(cstr != nullptr) {
      delay = atoi(cstr);
      std::cout << "READ_DELAY: " << delay << std::endl;
    } else delay = d;
  }
};

  //alignas(128) int timestamp_read::delay = 800;

struct alignas(64) timestamp_write {
  std::atomic<TS> stamp;
  alignas(128) static int delay;
  bool less(TS a, TS b) { return a < b;}
  bool equal(TS a, TS b) { return a == b;}
  
  TS get_stamp() {return stamp.load();}

  TS get_read_stamp() {return stamp.load();}

  TS get_write_stamp() {
    TS ts = stamp.load();

    if(delay == -1) {
      // delay to reduce contention
      for(volatile int i = 1; i < write_delay; i++) {}
      // std::atomic_thread_fence(std::memory_order_seq_cst);

      // only update timestamp if has not changed
      if (stamp.load() == ts) {
        if(stamp.fetch_add(1) == ts) {
          if(write_delay >= 2) write_delay /= 2;
        }
        else {
          if(write_delay < 256) write_delay *= 2;
        }
      }
    } else {
      for(volatile int i = 1; i < delay; i++) {}
      if (stamp.load() == ts)
        stamp.fetch_add(1);
    }
    return ts+1;
  }
  timestamp_write() : stamp(2) {
    auto cstr = std::getenv("WRITE_DELAY");
    if(cstr != nullptr)
      delay = atoi(cstr);
    else
      delay = -1;
    // std::cout << "WRITE_DELAY: " << delay << std::endl;
  }
};

alignas(128) int timestamp_write::delay = 200;

struct alignas(64) timestamp_multiple {
  static constexpr int slots = 4;
  static constexpr int gap = 16;
  static constexpr int delay = 300;
  std::atomic<TS> stamps[slots*gap];
  bool less(TS a, TS b) { return a < b;}
  bool equal(TS a, TS b) { return a == b;}
  
  inline TS get_write_stamp() {
    TS total = 0;
    for (int i = 0; i < slots; i++)
      total += stamps[i*gap].load();
    return total;
  }

  TS get_read_stamp() {
    TS ts = get_write_stamp();
    //int i = sched_getcpu() / 36; //% slots;
    int i = epoch::internal::worker_id() % slots;
    for(volatile int j = 0; j < delay ; j++);
    TS tsl = stamps[i*gap].load();
    if (ts == get_write_stamp()) {
      //int i = epoch::internal::worker_id() % slots;
      stamps[i*gap].compare_exchange_strong(tsl,tsl+1);
    }
    return ts;
  }

  timestamp_multiple() {
    for (int i = 0; i < slots; i++) stamps[i*gap] = 1;
  }
};


struct alignas(64) timestamp_no_inc {
  std::atomic<TS> stamp;
  bool less(TS a, TS b) { return a < b;}
  bool equal(TS a, TS b) { return a == b;}
  TS get_stamp() {return stamp.load();}
  TS get_read_stamp() {return stamp.load();}
  TS get_write_stamp() {return stamp.load();}
  timestamp_no_inc() : stamp(1) {}
};

// works well if mostly reads or writes
// if stamp is odd then in write mode, and if even in read mode
// if not in the right mode, then increment to put in the right mode
// thread_local float read_backoff = 50.0;
// thread_local float write_backoff = 1000.0;

thread_local int adaptive_backoff = 1;

struct alignas(64) timestamp_read_write {
  std::atomic<TS> stamp;
  alignas(128) static int delay;
  // unfortunately quite sensitive to the delays
  // these were picked empirically for a particular machine (aware)
  //static constexpr int write_delay = 1500;
  //static constexpr int read_delay = 150;

  TS get_stamp() {return stamp.load();}

  inline TS get_write_stamp() {
    TS s = stamp.load();
    if (s % 2 == 1) return s;
    if(delay == -1) {
      for(volatile int j = 1; j < adaptive_backoff ; j++) {}
    }
    else {
      for(volatile int i = 1; i < delay; i++) {}
    }
    if (s != stamp.load()) return s;
    TS tmp = s;
    if (stamp.compare_exchange_strong(tmp, s+1)) {
      if (adaptive_backoff >= 2) adaptive_backoff /= 2;
    } else if (adaptive_backoff < 256) adaptive_backoff *= 2;
    return s+1; // return new stamp
  }

  TS get_read_stamp() {
    TS s = stamp.load();
    if (s % 2 == 0) return s;
    if(delay == -1) {
      for(volatile int j = 1; j < adaptive_backoff ; j++) {}
    }
    else {
      for(volatile int j = 1; j < delay ; j++) {}
    }
    if (s != stamp.load()) return s;
    TS tmp = s;
    if (stamp.compare_exchange_strong(tmp, s+1)) {
      if (adaptive_backoff >= 2) adaptive_backoff /= 2;
    } else if (adaptive_backoff < 256) adaptive_backoff *= 2;
    return s; // return old stamp
  }

  TS current() { return stamp.load();}
  
  timestamp_read_write() : stamp(1) {
    auto cstr = std::getenv("WRITE_DELAY");
    // delay = 0;
    if(cstr != nullptr)
      delay = atoi(cstr);
    else
      delay = -1;
    // std::cout << "INC_DELAY: " << delay << std::endl;
  }
};

alignas(128) int timestamp_read_write::delay = 200;

// struct alignas(64) timestamp_read_write {
//   std::atomic<TS> stamp;
//   // unfortunately quite sensitive to the delays
//   // these were picked empirically for a particular machine (aware)
//   //static constexpr int write_delay = 1500;
//   //static constexpr int read_delay = 150;

//   TS get_stamp() {return stamp.load();}

//   inline TS get_write_stamp() {
//     TS s = stamp.load();
//     if (s % 2 == 1) return s;
//     for(volatile int j = 300; j < round(write_backoff) ; j++);
//     if (s != stamp.load()) return s;
//     TS tmp = s;
//     if (stamp.compare_exchange_strong(tmp, s+1)) {
//       if (write_backoff >= 100) write_backoff -= 200;
//     } else if (write_backoff < 600) write_backoff += 100;
//     return s+1; // return new stamp
//   }

//   TS get_read_stamp() {
//     TS s = stamp.load();
//     if (s % 2 == 0) return s;
//     for(volatile int j = 300; j < round(write_backoff) ; j++);
//     if (s != stamp.load()) return s;
//     TS tmp = s;
//     if (stamp.compare_exchange_strong(tmp, s+1)) {
//       if (write_backoff >= 100) write_backoff -= 200;
//     } else if (write_backoff < 600) write_backoff += 100;
//     return s; // return old stamp
//   }

//   TS current() { return stamp.load();}
  
//   timestamp_read_write() : stamp(1) {}
// };

#ifdef ReadStamp
timestamp_read global_stamp;
#elif TL2Stamp
  timestamp_tl2 global_stamp{100};
#elif LazyStamp
timestamp_read global_stamp{100};
#elif WriteStamp
timestamp_write global_stamp;
#elif HWStamp
timestamp_read_hw global_stamp;
#elif HWWriteStamp
timestamp_write_hw global_stamp;
#elif NoIncStamp
timestamp_no_inc global_stamp;
#else
timestamp_read_write global_stamp;
#endif

const TS zero_stamp = 1;
thread_local TS local_stamp{-1};

// this is updated by the epoch-based reclamation
// Whenever an epoch is incremented this is set to the stamp
// from the previous increment (which is now safe to collect).
TS done_stamp = global_stamp.get_stamp();
TS prev_stamp = global_stamp.get_stamp();
thread_local TS current_stamp;
  
  bool add_epoch_hooks() {
    epoch::internal::get_epoch().before_epoch_hooks.push_back([&] {
       current_stamp = global_stamp.get_stamp();});
    epoch::internal::get_epoch().after_epoch_hooks.push_back([&] {
	done_stamp = prev_stamp;
	prev_stamp = current_stamp;});
    return true;
  }

  bool junk = add_epoch_hooks();

template <typename F>
auto with_snapshot_internal(F f) {
  return flck::with_epoch([&] {
    local_stamp = global_stamp.get_read_stamp();
    if constexpr (std::is_void_v<std::invoke_result_t<F>>) {
      f();
      local_stamp = -1;
    } else {
      auto r = f();
      local_stamp = -1;
      return r;
    }
  });
}

thread_local bool aborted = false;

#ifndef LazyStamp
template <typename F>
auto with_snapshot(F f, bool unused_parameter=false) {
  return flck::with_epoch([&] {
    local_stamp = global_stamp.get_read_stamp();
    if constexpr (std::is_void_v<std::invoke_result_t<F>>) {
      f();
      local_stamp = -1;
    } else {
      auto r = f();
      local_stamp = -1;
      return r;
    }
  });
}
#else
thread_local bool speculative = false;
parlay::sequence<long> num_retries(epoch::internal::num_workers()*16, 0);
void print_retries() {
  std::cout << " retries = " << parlay::reduce(num_retries)
	    << ", final stamp = " << global_stamp.get_read_stamp() <<std::endl;
}

template <typename F>
auto with_snapshot(F f, bool use_speculative=true) {
  if(use_speculative) {
    return with_epoch([&] {
      local_stamp = global_stamp.get_stamp();
      aborted = false;
      speculative = true;
      if constexpr (std::is_void_v<std::invoke_result_t<F>>) {
        f();
        speculative = false;
        if (aborted) {
          aborted = false;
          num_retries[epoch::internal::worker_id()*16]++;
          global_stamp.increment_stamp(local_stamp);
          f();
        }
        local_stamp = -1;
      } else {
        auto r = f();
        speculative = false;
        if (aborted) {
          aborted = false;
          num_retries[epoch::internal::worker_id()*16]++;
          global_stamp.increment_stamp(local_stamp);
          r = f();
        }
        local_stamp = -1;
        return r;
      }
    });
  } else {
    return flck::with_epoch([&] {
      local_stamp = global_stamp.get_read_stamp();
      if constexpr (std::is_void_v<std::invoke_result_t<F>>) {
        f();
        local_stamp = -1;
      } else {
        auto r = f();
        local_stamp = -1;
        return r;
      }
    });
  }
}

#endif
} // namespace verlib




  
