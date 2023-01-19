#pragma once
#include <atomic>
#include "epoch.h"
#include "no_tagged.h"

namespace flck {

template <typename V>
struct atomic {
private:
  std::atomic<V> v;
public:
  static_assert(sizeof(V) <= 4 || std::is_pointer<V>::value,
    "Type for mutable must be a pointer or at most 4 bytes");
  atomic(V v) : v(v) {}
  atomic() : v(0) {}
  void init(V vv) {v = vv;}
  V load() {return v.load();}
  V read() {return v.load();}
  V read_snapshot() {return v.load();}
  V read_cur() {return v.load();}
  void store(V vv) { v = vv;}
  bool cas(V old_v, V new_v) {
    return (v.load() == old_v &&
	    v.compare_exchange_strong(old_v, new_v));}
  bool cas_ni(V old_v, V new_v) {
    return cas(old_v,new_v);}
  void cam(V old_v, V new_v) {
    if (v.load() == old_v)
      v.compare_exchange_strong(old_v,new_v);}
  V operator=(V b) {store(b); return b; }

  // compatibility with multiversioning
  void validate() {}
};

template <typename V>
using write_once = atomic<V>;

template <typename V>
struct atomic_write_once {
  std::atomic<V> v;
  atomic_write_once(V initial) : v(initial) {}
  atomic_write_once() {}
  V load() {return v.load();}
  V load_ni() {return v.load();}
  V read() {return v.load();}
  void init(V vv) { v = vv; }
  void store(V vv) { v = vv; }
  bool cas_ni(V exp_v, V new_v) {return v.compare_exchange_strong(exp_v, new_v);}
  V operator=(V b) { store(b); return b; }
  // inline operator V() { return load(); } // implicit conversion
};

template <typename T>
using memory_pool = internal::mem_pool<T>;

  // to make consistent with lock free implementation
  namespace internal {
    template <typename T>
    using tagged = no_tagged<T>;
  }

  template <typename F>
  bool skip_if_done(F f) {f(); return true;}

  template <typename F>
  bool skip_if_done_no_log(F f) {f(); return true;}

  template<typename V>
  V commit(V v) {return v;}

  template <typename F>
  void non_idempotent(F f) { f();}

} // namespace flck

