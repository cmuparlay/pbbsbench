#pragma once
#include "../flock/flock.h"
#include "timestamps.h"

namespace verlib {
  const TS tbd = std::numeric_limits<TS>::max()/4;

  
  template <typename T>
  using atomic = flck::atomic<T>;
  using lock = flck::lock;
  using atomic_bool = flck::atomic_write_once<bool>;
  using flck::memory_pool;
  template <typename F>
  auto do_now(F f) {return f();}

struct versioned {};

struct version_link {
  flck::atomic_write_once<long> time_stamp;
  version_link* next_version;
  void* value;
  version_link() : time_stamp(tbd) {}
  version_link(TS time, version_link* next, void* value) :
    time_stamp(time), next_version(next), value(value) {}  
};

flck::memory_pool<version_link> link_pool;

template <typename V>
struct versioned_ptr {
private:
  flck::atomic<version_link*> v;

  static version_link* set_stamp(version_link* ptr) {
    if (ptr->time_stamp.load_ni() == tbd) {
      TS old_t = tbd;
      TS new_t = global_stamp.get_write_stamp();
      if (ptr->time_stamp.load_ni() == tbd)
        ptr->time_stamp.cas_ni(old_t, new_t);
    }
    return ptr;
  }

  static version_link* init_ptr(V* ptr) {
    return link_pool.New(zero_stamp, nullptr, (void*) ptr);
  }

  bool idempotent_cas(version_link* old_v, version_link* new_v) {
#ifdef NoHelp
    return v.cas(old_v, new_v);
#else
    v.cam(old_v, new_v);
    return (v.load() == new_v ||
	    new_v->time_stamp.load() != tbd);
#endif
  }
  
public:

  versioned_ptr(): v(init_ptr(nullptr)) {}
  versioned_ptr(V* ptr) : v(init_ptr(ptr)) {}
  ~versioned_ptr() { link_pool.Delete(v.load()); }
  void init(V* ptr) {v = init_ptr(ptr);}
  
  V* read_snapshot() {
    version_link* head = set_stamp(v.load());
    while (global_stamp.less(local_stamp, head->time_stamp.load()))
      head = head->next_version;
#ifdef LazyStamp
    if (global_stamp.equal(head->time_stamp.load(), local_stamp) && speculative)
      aborted = true;
#endif
    return (V*) head->value;
  }

  V* load() {  // can be used anywhere
    if (local_stamp != -1) return read_snapshot();
    else return (V*) set_stamp(v.load())->value;
  }

  // only safe on journey
  V* read() {  return (V*) v.read()->value; }

  void validate() { set_stamp(v.load()); }

  void store(V* ptr) {
    version_link* old_link = v.load();
    version_link* new_link = link_pool.New(tbd, old_link, (void*) ptr);
    v = new_link;
    set_stamp(new_link);
    link_pool.Retire(old_link);    
  }

  bool cas(V* old_v, V* new_v) {
    version_link* old_link = v.load();
    if(old_link != nullptr) set_stamp(old_link);
    if (old_v != old_link->value) return false;
    if (old_v == new_v) return true;
    version_link* new_link = link_pool.New(tbd, old_link, (void*) new_v);
    if (idempotent_cas(old_link, new_link)) {
      set_stamp(new_link);
      link_pool.Retire(old_link);
      return true;
    }
    set_stamp(v.load());
    link_pool.Delete(new_link);
    return false;
  }

  V* operator=(V* b) {store(b); return b; }
};

  template <typename T, typename E>
  bool validate(flck::lock& lck, T* v, E expected) {
    return true;
  }

} // namespace verlib
