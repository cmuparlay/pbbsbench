#include <atomic>
#include <vector>
#include <limits>

#include "parlay/alloc.h"
#include "parlay/primitives.h"

#ifndef PARLAY_EPOCH_H_
#define PARLAY_EPOCH_H_

#ifndef NDEBUG
// Checks for corruption of bytes before and after allocated structures, as well as double frees.
// Requires some extra memory to pad the front and back of a structure.
#define EpochMemCheck 1
#endif

// Can use malloc instead of parlay::type_allocator
//#define USE_MALLOC 1

// If defined allows with_epoch to be nested (only outermost will set the epoch).
// Incurs slight overhead due to extra test, but allows wrapping a with_epoch
// around multiple operations which each do a with_epoch.
// This saves time since setting epoch number only needs to be fenced once on the outside.
//#define NestedEpochs 1

// Supports before_epoch_hooks and after_epoch_hooks, which are thunks
// that get run just before incrementing the epoch number and just after.
// The user set them with:
//    flck::internal::epoch.before_epoch_hooks.push_back(<mythunk>);
//    flck::internal::epoch.after_epoch_hooks.push_back(<myotherthunk>);


// ***************************
// epoch structure
// ***************************

namespace epoch {
namespace internal {
  // thread_local int __count = 0;

  inline int worker_id() {
    // if(__count++ % 1024 == 0)
      // std::cerr << "my_thread_id: " << parlay::my_thread_id() << std::endl;
    return parlay::my_thread_id();
  }

  inline int num_workers() {
    if (const auto env_p = std::getenv("PARLAY_NUM_THREADS")) {
      return std::stoi(env_p);
    } else {
      return std::thread::hardware_concurrency()+2;
    }
    // std::cerr << "num_thread_ids: " << parlay::num_thread_ids() << std::endl;
    // if(parlay::num_thread_ids() == 0) abort();
    // return parlay::num_thread_ids();
  }

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
    int workers = num_workers();
    announcements = std::vector<announce_slot>(workers);
    current_epoch = 0;
  }

  void print_announce() {
    for(auto& ann : announcements)
      std::cout << ann.last << " ";
    std::cout << std::endl;
  }

  void clear_announce() {
    for(auto& ann : announcements)
      ann.last = -1;
  }

  long get_current() {
    return current_epoch.load();
  }
  
  long get_my_epoch() {
    return announcements[worker_id()].last;
  }

  void set_my_epoch(long e) {
    announcements[worker_id()].last = e;
  }

  std::pair<bool,int> announce() {
    size_t id = worker_id();
    while (true) {
      long current_e = get_current();
#ifdef NestedEpochs
      long old = -1l;
      bool succeeded = false;
      if(announcements[id].last.load() == old)
        succeeded = announcements[id].last.compare_exchange_strong(old, current_e);
      // if(!succeeded) abort();
      if (get_current() == current_e) return std::pair(succeeded, id);
#else
      long tmp = current_e;
      // apparently an exchange is faster than a store (write and fence)
      announcements[id].last.exchange(tmp, std::memory_order_seq_cst);
      if (get_current() == current_e) return std::pair(true, id);
#endif
    }
  }

  void unannounce(size_t id) {
    assert(announcements[id].last.load() != -1l);
    announcements[id].last.store(-1l, std::memory_order_release);
  }

  void update_epoch() {
    size_t id = worker_id();
    int workers = num_workers();
    long current_e = get_current();
    bool all_there = true;
    // check if everyone is done with earlier epochs
    for (int i=0; i < workers; i++)
      if ((announcements[i].last != -1l) && announcements[i].last < current_e) {
        all_there = false;
        break;
      }
    // if so then increment current epoch
    if (all_there) {
      for (auto h : before_epoch_hooks) h();
      if (current_epoch.compare_exchange_strong(current_e, current_e+1)) {
        for (auto h : after_epoch_hooks) h();
      }
    }
  }

};

  extern inline epoch_s& get_epoch() {
    static epoch_s epoch;
    return epoch;
  }

// ***************************
// epoch pools
// ***************************

struct Link {
  Link* next;
  bool skip;
  void* value;
};

  // x should point to the skip field of a link
  inline void undo_Retire(bool* x) { *x = true;}
  inline void undo_allocate(bool* x) { *x = false;}

#ifdef USE_MALLOC
  inline Link* allocate_link() {return (Link*) malloc(sizeof(Link));}
  inline void free_link(Link* x) {return free(x);}
#else
  using list_allocator = typename parlay::type_allocator<Link>;
  inline Link* allocate_link() {return list_allocator::alloc();}
  inline void free_link(Link* x) {return list_allocator::free(x);}
#endif
  
  using namespace std::chrono;

template <typename xT>
struct alignas(64) memory_pool {
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
    auto i = worker_id();
    auto &pid = pools[i];
    advance_epoch(i, pid);
    Link* lnk = allocate_link();
    lnk->next = pid.current;
    lnk->value = p;
    lnk->skip = false;
    pid.current = lnk;
    return &(lnk->skip);
  }

  // destructs and frees a linked list of objects 
  void clear_list(Link* ptr) {
    // abort();
    while (ptr != nullptr) {
      Link* tmp = ptr;
      ptr = ptr->next;
      if (!tmp->skip) {
#ifdef EpochMemCheck
        paddedT* x = pad_from_T((T*) tmp->value);
        if (x->head != 10 || x->tail != 10) {
          if (x->head == 55) std::cerr << "double free" << std::endl;
          else std::cerr << "corrupted head" << std::endl;
          if (x->tail != 10) std::cerr << "corrupted tail" << std::endl;
          assert(false);
        }
#endif
        Delete((T*) tmp->value);
      }
      free_link(tmp);
    }
  }

  // computes size of list
  long size_of(Link* ptr) {
    long sum = 0;
    while (ptr != nullptr) {
      Link* tmp = ptr;
      ptr = ptr->next;
      sum++;
    }
    return sum;
  }

  void advance_epoch(int i, old_current& pid) {
    epoch_s& epoch = get_epoch();
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

#ifdef  EpochMemCheck
  using nodeT = paddedT;
#else
  using nodeT = xT;
#endif

#ifdef USE_MALLOC
  nodeT* allocate_node() {return (nodeT*) malloc(sizeof(nodeT));}
  void free_node(nodeT* x) {return free(x);}
#else
  using Allocator = parlay::type_allocator<nodeT>;
  nodeT* allocate_node() { return Allocator::alloc();}
  void free_node(nodeT* x) { return Allocator::free(x);}
#endif
  
public:
  using T = xT;
  
  memory_pool() {
    workers = num_workers();
    update_threshold = 10 * workers;
    pools = std::vector<old_current>(workers);
    for (int i = 0; i < workers; i++) {
      pools[i].count = parlay::hash64(i) % update_threshold;
      pools[i].time = system_clock::now();
    }
  }

  memory_pool(const memory_pool&) = delete;
  ~memory_pool() {} // clear(); }

  // noop since epoch announce is used for the whole operation
  void acquire(T* p) { }
  
  paddedT* pad_from_T(T* p) {
     size_t offset = ((char*) &((paddedT*) p)->value) - ((char*) p);
     return (paddedT*) (((char*) p) - offset);
  }
  
  // destructs and frees the object immediately
  void Delete(T* p) {
     p->~T();
#ifdef EpochMemCheck
     paddedT* x = pad_from_T(p);
     x->head = 55;
     free_node(x);
#else
     free_node(p);
#endif
  }

  template <typename ... Args>
  T* New(Args... args) {
#ifdef EpochMemCheck
    paddedT* x = allocate_node();
    x->pad = x->head = x->tail = 10;
    T* newv = &x->value;
    new (newv) T(args...);
    assert(check_not_corrupted(newv));
#else
    T* newv = allocate_node();
    new (newv) T(args...);
#endif
    return newv;
  }

  bool check_not_corrupted(T* ptr) {
#ifdef EpochMemCheck
    paddedT* x = pad_from_T(ptr);
    if (x->pad != 10) std::cerr << "memory_pool: pad word corrupted" << std::endl;
    if (x->head != 10) std::cerr << "memory_pool: head word corrupted" << std::endl;
     if (x->tail != 10) std::cerr << "memory_pool: tail word corrupted" << std::endl;
    return (x->pad == 10 && x->head == 10 && x->tail == 10);
#endif
    return true;
  }
                           
  template <typename F, typename ... Args>
  // f is a function that initializes a new object before it is shared
  T* new_init(const F& f, Args... args) {
    T* x = New(args...);
    f(x);
    return x;
  }

  // retire and return a pointer if want to undo the retire
  bool* Retire(T* p) {
    return add_to_current_list((void*) p);}
  
  // clears all the lists 
  // to be used on termination
  void clear() {
    get_epoch().update_epoch();
    for (int i=0; i < pools.size(); i++) {
      clear_list(pools[i].old);
      clear_list(pools[i].current);
      pools[i].old = pools[i].current = nullptr;
    }
  }

  void reserve(size_t n) {
#ifndef USE_MALLOC
    Allocator::reserve(n);
#endif
  }

  void stats() {
    get_epoch().print_announce();
    // get_epoch().clear_announce();
    std::cout << "epoch number: " << get_epoch().get_current() << std::endl;
    for (int i=0; i < pools.size(); i++) {
      std::cout << "pool[" << i << "] = " << size_of(pools[i].old) << ", " << size_of(pools[i].current) << std::endl;
    }
#ifndef USE_MALLOC
    Allocator::print_stats();
#endif
  }

  void shuffle(size_t n) {}
    
};

template <typename T>
extern inline memory_pool<T>& get_pool() {
  static memory_pool<T> pool;
  return pool;
}

} // end namespace internal

template <typename Thunk>
auto with_epoch(Thunk f) {
  auto& epoch = internal::get_epoch();
  auto [not_in_epoch, id] = epoch.announce();
  if constexpr (std::is_void_v<std::invoke_result_t<Thunk>>) {
    f();
#ifdef NestedEpochs
    if (not_in_epoch) 
#endif
      epoch.unannounce(id);
  } else {
    auto v = f();
#ifdef NestedEpochs
    if (not_in_epoch) 
#endif
      epoch.unannounce(id);
    return v;
  }
}

  template <typename T>
  using memory_pool = internal::memory_pool<T>;

  template <typename T>
  struct memory_pool_ {
    template <typename ... Args>
    static T* New(Args... args) {
      return internal::get_pool<T>().New(std::forward<Args>(args)...);}
    static void Delete(T* p) {internal::get_pool<T>().Delete(p);}
    static bool* Retire(T* p) {return internal::get_pool<T>().Retire(p);}
    static void clear() {internal::get_pool<T>().clear();}
  };

} // end namespace flck

#endif //PARLAY_EPOCH_H_
