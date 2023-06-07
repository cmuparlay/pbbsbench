#include <atomic>
#include "tagged.h"
#include "lf_log.h"

#pragma once

namespace flck {

template <typename V>
struct atomic {
private:
  using IT = size_t;
  using TV = internal::tagged<V>;

  IT get_val(internal::Log &p) {
    return p.commit_value(v.load()).first; }

public:
  std::atomic<IT> v;
  static_assert(sizeof(V) <= 4 || std::is_pointer<V>::value,
    "Type for mutable must be a pointer or at most 4 bytes");

  // not much to it.  heavy lifting done in TV
  atomic(V vv) : v(TV::init(vv)) {}
  atomic() : v(TV::init(0)) {}
  void init(V vv) {v = TV::init(vv);}
  V load() {return TV::value(get_val(internal::lg));}
  V load_ni() {return TV::value(v.load());}
  V read() {return TV::value(v.load());}
  V read_snapshot() {return TV::value(v.load());}
  void store(V vv) {TV::cas(v, get_val(internal::lg), vv);}
  bool cas(V old_v, V new_v) { // not safe inside locks
    assert(internal::lg.is_empty());
    return cas_ni(old_v, new_v);
  }
  bool cas_ni(V old_v, V new_v) { 
    IT old_t = v.load();
    return (TV::value(old_t) == old_v &&
      TV::cas(v, old_t, new_v, true));}
  void cam(V oldv, V newv) {
    IT old_t = get_val(internal::lg);
    if (TV::value(old_t) == oldv)
      TV::cas(v, old_t, newv);}
  V operator=(V b) {store(b); return b; }

  // compatibility with multiversioning
  void validate() {}
  
  // operator V() { return load(); } // implicit conversion
};

template <typename V>
struct atomic_double {
private:
  struct alignas(16) TV {size_t count; V val; };
  TV v;
  void cam(TV oldv, TV newv) {
    __int128 oldvi = *((__int128*) &oldv);
    __int128 newvi = *((__int128*) &newv);
    __sync_bool_compare_and_swap((__int128*) this, oldvi, newvi);
  }
      
public:
  static_assert(sizeof(V) <= 8,
    "Type for mutable_double must be at most 8 bytes");

  atomic_double(V vv) : v({1, vv}) {}
  atomic_double() : v({1, 0}) {}
  V load() {return internal::lg.commit_value_safe(v.val).first;}
  V read() {return v.val;}
  void init(V vv) {v.val = vv;}
  void store(V newv) {
    size_t cnt = internal::lg.commit_value(v.count).first;
#ifdef NoSkip  // used to test performance without optimization
    cam({cnt, v.val}, {cnt+1, newv});
#else
    // skip if done for efficiency
    skip_if_done_no_log([&] { cam({cnt, v.val}, {cnt+1, newv});});
#endif
  }
  V operator=(V b) {
    store(b);
    return b;
  }
};

template <typename V>
struct atomic_aba_free {
private:
   constexpr static unsigned long set_bit= (1ul << 63);
public:
  std::atomic<V> v;
  atomic_aba_free(V initial) : v(initial) {}
  atomic_aba_free() {}
  V load() { // set then mask high bit to ensure not zero
    size_t x = internal::lg.commit_value((size_t) v.load() | set_bit).first;
    return (V) (x & ~set_bit);
  }
  void init(V vv) { v = vv; }
  void store(V vv) {
    V old_v = load();
    v.compare_exchange_strong(old_v, vv);
  }
  void cam(V expected, V new_v) {
    V old_v = load();
    if (expected == old_v)
      v.compare_exchange_strong(old_v, new_v);
  }
  V operator=(V b) { store(b); return b; }
};

template <typename V>
struct atomic_write_once {
private:
   constexpr static unsigned long set_bit= (1ul << 63);
public:
  std::atomic<V> v;
  atomic_write_once(V initial) : v(initial) {}
  atomic_write_once() {}
  V load() { // set then mask high bit to ensure not zero
    size_t x = internal::lg.commit_value((size_t) v.load() | set_bit).first;
    return (V) (x & ~set_bit);
  }
  V load_ni() {return v.load();}
  void init(V vv) { v = vv; }
  void store(V vv) { v = vv; }
  bool cas_ni(V exp_v, V new_v) {return v.compare_exchange_strong(exp_v, new_v);}
  V operator=(V b) { store(b); return b; }
  // inline operator V() { return load(); } // implicit conversion
};

// *****************************
// Memory pool using epoch based collection and safe (idempotent)
// allocation and retire in a lock.
// *****************************

 namespace internal {
   struct lock;
 }
 
template <typename T, typename Pool=internal::mem_pool<T>>
struct memory_pool {
  Pool pool;

  void reserve(size_t n) { pool.reserve(n);}
  void clear() { pool.clear(); }
  void stats() { pool.stats();}
  void shuffle(size_t n) { pool.shuffle(n);}
  
  void acquire(T* p) { pool.acquire(p);}
  
  bool* retire(T* p) {
    assert(p != nullptr);
    if (!internal::helping)
      return internal::with_empty_log([&] {return pool.retire(p);});
    else return nullptr;
    //auto x = internal::lg.commit_value_safe(p);
    //if (x.second) // only retire if first try
    //  internal::with_empty_log([&] {pool.retire(p);}); 
  }

  bool* retire_ni(T* p) {
    assert(p != nullptr);
    return internal::with_empty_log([&] {return pool.retire(p);});
  }
	      
  void destruct(T* p) {
    assert(p != nullptr);
    auto x = internal::lg.commit_value_safe(p);
    if (x.second) // only retire if first try
      internal::with_empty_log([&] {pool.destruct(p);}); 
  }

  void destruct_ni(T* p) {
    assert(p != nullptr);
    internal::with_empty_log([&] {pool.destruct(p);}); 
  }

  template <typename F, typename ... Args>
  // f is a function that initializes a new object before it is shared
  T* new_init(F f, Args... args) {
    //run f without logging (i.e. an empty log)
    T* newv = internal::with_log(internal::Log(), [&] { 
	T* x = pool.new_obj(args...);
	f(x);
	return x;});
    auto r = internal::lg.commit_value(newv);
    if (!r.second)  // destruct if already initialized
      internal::with_empty_log([&] {pool.destruct(newv);});
    return r.first;
  }

  // Idempotent allocation
  template <typename ... Args>
  T* new_obj(Args ...args) {
    return new_obj_fl(args...).first;
  }

protected:
  friend class internal::lock;

  // The following protected routines are only used internally
  // in the lock code (not accessible to the user)
  
  // Returns a pointer to the new object and a possible pointer
  // to a location in the log containing the pointer.
  // The location is null if this was not the first among thunks
  // to allocate the object.
  // The returned pointer can be one of done_true or done_false
  // if the object is already retired using retire_acquired.
  template <typename ... Args>
  std::pair<T*,internal::log_entry*> new_obj_acquired(Args... args) {
    auto [ptr,fl] = new_obj_fl(args...);
    if (internal::lg.is_empty()) return std::pair(ptr, nullptr);
    internal::log_entry* l = internal::lg.current_entry();
    if (!fl && !is_done(ptr)) {
      pool.acquire(ptr);
      return std::make_pair((T*) l->load(), nullptr);
    } else return std::make_pair(ptr, fl ? l : nullptr);
  }

  // le must be a value returned as the second return value of new_obj_acquired.
  // It will be either be null or a pointer to a log entry containing p.
  // If non-null then it clears p from the log by replacing it with
  // the result (true or false) so that p can be safely reclaimed.
  // It then retires p.
  // It is important that only one of the helping thunks is passed an le that
  // is not null, otherwise could be retired multiple times
  template<typename TT>
  void retire_acquired_result(T* p, internal::log_entry* le, std::optional<TT> result) {
    if (internal::lg.is_empty()) pool.retire(p);
    else if (le != nullptr) {
      *le = tag_result(result);
      internal::with_empty_log([&] { pool.retire(p);});
    }
  }

  bool is_done(T* p) {return is_done_flag(p);}

  template <typename RT>
  std::optional<RT> done_val_result(T* p) {
    auto r = extract_result(p);
    if (r.has_value()) return (RT) r.value();
    else return {};
  }
  
private:
  bool done_val(T* p) {return extract_bool(p);}

  // this version also returns a flag to specify whether actually allocated
  // vs., having been allocated by another instance of a thunk
  template <typename ... Args>
  std::pair<T*,bool> new_obj_fl(Args... args) {
    // TODO: helpers might do lots of allocates and frees,
    // can potentially optimize by checking if a log value has already been committed.
    T* newv = internal::with_log(internal::Log(),
				 [&] {return pool.new_obj(args...);});
    auto r = internal::lg.commit_value(newv);
    // if already allocated return back to pool
    if (!r.second) {
      internal::with_empty_log([&] {pool.destruct(newv);}); }
    return r;
  }

  // the following tags a long entry with the return value of a thunk
  // 1 at the 48th bit is true, 2 at the 48th bit is false
  
  bool is_done_flag(T* p) {
    return (((size_t) p) >> 48) > 0;
  }

  void* tag_bool(bool result) {
    return (void*) (result ? (1ul << 48) : (2ul << 48));
  }

  bool extract_bool(T* p) {
    return (((size_t) p) >> 48) == 1ul;
  }

  // a poor mans "optional".  The flag at the 48th bit indicates presence,
  // and the lower 48 bits are the value if present.
  template<typename TT>
  void* tag_result(std::optional<TT> result) {
    if(!result.has_value()) return (void*) (2ul << 48);
    else return (void*) ((1ul << 48) | (size_t) result.value());
  }

  std::optional<size_t> extract_result(T* p) {
    if (extract_bool(p))
      return ((size_t) p) & ((1ul << 48) - 1);
    return {};
  }
  
};

template<typename V>
V commit(V v) {return internal::lg.commit_value(v).first;}

template <typename F>
bool skip_if_done(F f) { return internal::skip_if_done(f); }

template <typename F>
bool skip_if_done_no_log(F f) { return internal::skip_if_done_no_log(f); }

template <typename F>
void non_idempotent(F f) { internal::with_empty_log(f); }

} // namespace flck

