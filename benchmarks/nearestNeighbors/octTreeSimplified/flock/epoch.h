#include <atomic>
#include <vector>
#include <limits>
#include <parlay/alloc.h>
#include <parlay/random.h>
#include <parlay/primitives.h>
#include <unordered_set>
//#include "timestamps.h"

#pragma once

#ifndef NDEBUG
  #define EpochMemCheck 1
#endif

// ***************************
// epoch structure
// ***************************

namespace flck {
namespace internal {

struct alignas(64) epoch_s {
	
  // functions to run when epoch is incremented
  std::vector<std::function<void()>> before_epoch_hooks;
  std::vector<std::function<void()>> after_epoch_hooks;
  
  struct alignas(64) announce_slot {
    std::atomic<long> last;
    announce_slot() : last(-1l) {}
  };

  std::vector<announce_slot> announcements;
  std::atomic<long> current_epoch;
  epoch_s() {
    int workers = parlay::num_workers();
    announcements = std::vector<announce_slot>(workers);
    current_epoch = 0;
  }

  long get_current() {
    return current_epoch.load();
  }

  long get_my_epoch() {
    size_t id = parlay::worker_id();
    return announcements[id].last;
  }

  void set_my_epoch(long e) {
    size_t id = parlay::worker_id();
    announcements[id].last = e;
  }

  void announce() {
    size_t id = parlay::worker_id();
    long current_e = get_current();
    // apparently an exchange is faster than a store (write and fence)
    announcements[id].last.exchange(current_e, std::memory_order_acquire);
  }

  void unannounce() {
    size_t id = parlay::worker_id();
    announcements[id].last.store(-1l, std::memory_order_release);
  }

  void update_epoch() {
    size_t id = parlay::worker_id();
    int workers = parlay::num_workers();
    long current_e = get_current();
    bool all_there = true;
    // check if everyone is done with earlier epochs
    for (int j=0; j<2; j++) //do twice
      for (int i=0; i < workers; i++)
      	if ((announcements[i].last != -1l) && announcements[i].last < current_e) 
      	  all_there = false;
    // if so then increment current epoch
    if (all_there) {
      for (auto h : before_epoch_hooks) h();
      if (current_epoch.compare_exchange_strong(current_e, current_e+1)) {
	for (auto h : after_epoch_hooks) h();
      }
    }
  }
};

epoch_s epoch;

// ***************************
// epoch pools
// ***************************

struct Link {
  Link* next;
  bool skip;
  void* value;
};

  // x should point to the skip field of a link
  void undo_retire(bool* x) { *x = true;}
  void undo_allocate(bool* x) { *x = false;}
  
using list_allocator = parlay::type_allocator<Link>;

  using namespace std::chrono;

template <typename xT>
struct alignas(64) mem_pool {
private:

  static constexpr double milliseconds_between_epoch_updates = 20.0;
  long update_threshold;
  using sys_time = time_point<std::chrono::system_clock>;

  // each thread keeps one of these
  struct alignas(256) old_current {
    Link* old;  // linked list of retired items from previous epoch
    Link* current; // linked list of retired items from current epoch
    long epoch; // epoch on last retire, updated on a retire
    long count; // number of retires so far, reset on updating the epoch
    sys_time time; // time of last epoch update
    old_current() : old(nullptr), current(nullptr), epoch(0) {}
  };

  // only used for debugging (i.e. EpochMemCheck=1).
  struct paddedT {
    long pad;
    std::atomic<long> head;
    xT value;
    std::atomic<long> tail;
  };

  std::vector<old_current> pools;
  int workers;

  bool* add_to_current_list(void* p) {
    auto i = parlay::worker_id();
    auto &pid = pools[i];
    advance_epoch(i, pid);
    Link* lnk = list_allocator::alloc();
    lnk->next = pid.current;
    lnk->value = p;
    lnk->skip = false;
    pid.current = lnk;
    return &(lnk->skip);
  }

  // destructs and frees a linked list of objects 
  void clear_list(Link* ptr) {
    while (ptr != nullptr) {
      Link* tmp = ptr;
      ptr = ptr->next;
      if (!tmp->skip) {
#ifdef EpochMemCheck
	paddedT* x = pad_from_T((T*) tmp->value);
	if (x->head != 10 || x->tail != 10) {
	  if (x->head == 55) std::cout << "double free" << std::endl;
	  else std::cout << "corrupted head" << std::endl;
	  if (x->tail != 10) std::cout << "corrupted tail" << std::endl;
	  assert(false);
	}
#endif
	destruct((T*) tmp->value);
      }
      list_allocator::free(tmp);
    }
  }

  void advance_epoch(int i, old_current& pid) {
    if (pid.epoch + 1 < epoch.get_current()) {
      clear_list(pid.old);
      pid.old = pid.current;
      pid.current = nullptr;
      pid.epoch = epoch.get_current();
    }
    // a heuristic
    auto now = system_clock::now();
    if (++pid.count == update_threshold  ||
	duration_cast<milliseconds>(now - pid.time).count() >
	milliseconds_between_epoch_updates * (1 + ((float) i)/workers)) {
      pid.count = 0;
      pid.time = now;
      epoch.update_epoch();
    }
  }


public:
  using T = xT;
#ifdef  EpochMemCheck
  using Allocator = parlay::type_allocator<paddedT>;
#else
  using Allocator = parlay::type_allocator<T>;
#endif
  
  mem_pool() {
    workers = parlay::num_workers();
    update_threshold = 10 * workers;
    pools = std::vector<old_current>(workers);
    for (int i = 0; i < workers; i++) {
      pools[i].count = parlay::hash64(i) % update_threshold;
      pools[i].time = system_clock::now();
    }
  }

  mem_pool(const mem_pool&) = delete;
  ~mem_pool() { clear(); }

  // noop since epoch announce is used for the whole operation
  void acquire(T* p) { }
  
  void reserve(size_t n) { Allocator::reserve(n);}
  void stats() { Allocator::print_stats();}

  void shuffle(size_t n) {
    n = std::max(n, 1000000ul);
    auto ptrs = parlay::tabulate(n, [&] (size_t i) {return Allocator::alloc();});
    ptrs = parlay::random_shuffle(ptrs);
    parlay::parallel_for(0, n, [&] (size_t i) {Allocator::free(ptrs[i]);});
  }

  paddedT* pad_from_T(T* p) {
     size_t offset = ((char*) &((paddedT*) p)->value) - ((char*) p);
     return (paddedT*) (((char*) p) - offset);
  }
  
  // destructs and frees the object immediately
  void destruct(T* p) {
     p->~T();
#ifdef EpochMemCheck
     paddedT* x = pad_from_T(p);
     x->head = 55;
     Allocator::free(x);
#else
     Allocator::free(p);
#endif
  }

  template <typename ... Args>
  T* new_obj(Args... args) {
#ifdef EpochMemCheck
    paddedT* x = Allocator::alloc();
    x->head = x->tail = 10;
    T* newv = &x->value;
#else
    T* newv = Allocator::alloc();
#endif
    new (newv) T(args...);
    return newv;
  }

  template <typename F, typename ... Args>
  // f is a function that initializes a new object before it is shared
  T* new_init(F f, Args... args) {
    T* x = new_obj(args...);
    f(x);
    return x;
  }

  // retire and return a pointer if want to undo the retire
  bool* retire(T* p) {
    return add_to_current_list((void*) p);}
  
  // clears all the lists and terminates the underlying allocator
  // to be used on termination
  void clear() {
    epoch.update_epoch();
    for (int i=0; i < pools.size(); i++) {
      clear_list(pools[i].old);
      clear_list(pools[i].current);
      pools[i].old = pools[i].current = nullptr;
    }
    Allocator::finish();
  }
};

} // end namespace internal

template <typename Thunk>
auto with_epoch(Thunk f) {
  internal::epoch.announce();
  if constexpr (std::is_void_v<std::invoke_result_t<Thunk>>) {
    f();
    internal::epoch.unannounce();
  } else {
    auto v = f();
    internal::epoch.unannounce();
    return v;
  }
}
} // end namespace flck
