#pragma once

#include "flock/flock.h"
#include "timestamps.h"

namespace vl {
  
using IT = size_t;

parlay::sequence<long> i_counts(parlay::num_workers()*16, 0);
void print_counts() {
  std::cout << " indirect = " << parlay::reduce(i_counts) << std::endl;
}

struct versioned {
  size_t foo;
  std::atomic<TS> time_stamp;
  versioned* next_version;
  static constexpr size_t init_ptr =(1ul << 48) - 2;
  versioned* add_tag(versioned* v, bool tag) {
    return (versioned*) ((size_t) v + tag);}
  bool is_indirect() {return (IT) next_version & 1;}
  versioned* get_next() {return (versioned*) ((IT) next_version & ~1ul);}
  TS read_stamp() {return time_stamp.load();}
  TS load_stamp() {return flck::commit(time_stamp.load());}
  void set_stamp(TS t) {
    TS old = tbd;
    if(time_stamp.load() == tbd)
      time_stamp.compare_exchange_strong(old, t);
  }
  versioned() : time_stamp(tbd), next_version((versioned*) init_ptr) {}
  versioned(versioned* next, bool is_indirect)
    : time_stamp(tbd), next_version(add_tag(next,is_indirect)) {}
};

struct plink : versioned {
  void* value;
  plink(versioned* next, void* value) : versioned{next, true}, value(value) {}  
};

 flck::memory_pool<plink> link_pool;

template <typename V>
struct versioned_ptr {
private:
  flck::atomic<V*> v;

  static V* set_stamp(V* ptr) {
    if (ptr != nullptr && ptr->read_stamp() == tbd)
      ptr->set_stamp(global_stamp.get_write_stamp());
    return ptr;
  }

  static V* set_zero_stamp(V* ptr) {
    if (ptr != nullptr && ptr->read_stamp() == tbd)
      ptr->time_stamp = zero_stamp;
    return ptr;
  }

  void shortcut(plink* ptr) {
    if (ptr->read_stamp() <= done_stamp)
      flck::non_idempotent([&] {
	  if (v.cas((V*) ptr, (V*) ptr->value))
	    link_pool.retire(ptr);});      
  }

  V* get_ptr(V* ptr) {
    if (ptr != nullptr && ptr->is_indirect()) {
#ifndef NoShortcut
      shortcut((plink*) ptr);
#endif
      return (V*) ((plink*) ptr)->value;
    } else return ptr;
  }

public:

  versioned_ptr(): v(0) {}
  versioned_ptr(V* ptr) : v(set_zero_stamp(ptr)) {}

  ~versioned_ptr() {
    plink* ptr = (plink*) v.load();
    if (ptr != nullptr && ptr->is_indirect())
      link_pool.destruct(ptr);
  }

  void init(V* ptr) {v = set_zero_stamp(ptr);}

  V* read_snapshot() {
    TS ls = local_stamp;
    V* head = v.read();
    set_stamp(head);
    while (head != nullptr && head->read_stamp() > ls)
      head = (V*) head->get_next();
#ifdef LazyStamp
    if (head != nullptr && head->read_stamp() == ls && speculative)
      aborted = true;
#endif
    return (((head != nullptr) && head->is_indirect()) ?
	    (V*) ((plink*) head)->value :
	    head);
  }

  V* load() {  // can be used anywhere
    if (local_stamp != -1) return read_snapshot();
    else return get_ptr(set_stamp(v.load()));
  }
  
  V* read() {  // only safe on journey
    return get_ptr(v.read());
  }

  void validate() {
    set_stamp(v.load());     // ensure time stamp is set
  }
  
  void store(V* ptr) {
    V* old_v = v.load();
    V* new_v = ptr;
    bool use_indirect = (ptr == nullptr || ptr->load_stamp() != tbd);

    if (use_indirect)
      new_v = (V*) link_pool.new_obj((versioned*) old_v, ptr);
    else ptr->next_version = old_v;

#ifdef NoShortcut
    v = new_v;
    if (old_v != nullptr && old_v->is_indirect()) 
      link_pool.retire((plink*) old_v);
#else
    v.cam(old_v, new_v);
    if (old_v != nullptr && old_v->is_indirect()) {
      // link_pool.retire((plink*) old_v);
      V* val = v.load();
      if (val != (V*) ((plink*) old_v)->value)
	      link_pool.retire((plink*) old_v);
      else v.cam(val, new_v);
    }
#endif
    set_stamp(new_v);
#ifndef NoShortcut
    if (use_indirect) shortcut((plink*) new_v);
#endif
  }
  
  bool cas(V* expv, V* newv) {
#ifndef NoShortcut
    for(int ii = 0; ii < 2; ii++) {
#endif
      V* new_v = newv;
      V* oldv = v.load();
      if(oldv != nullptr) set_stamp(oldv);
      if(get_ptr(oldv) != expv) return false;
      if(oldv == newv) return true;
      bool use_indirect = (newv == nullptr || newv->load_stamp() != tbd);

      if(use_indirect)
        new_v = (V*) link_pool.new_obj((versioned*) oldv, newv);
      else newv->next_version = oldv;

      bool succeeded = v.cas_ni(oldv, new_v);

      if(succeeded) {
        set_stamp(new_v);
        if(oldv != nullptr && oldv->is_indirect()) 
          link_pool.retire((plink*) oldv);
#ifndef NoShortcut
          if (use_indirect) shortcut((plink*) new_v);
#endif
        return true;
      }
      if(use_indirect) link_pool.destruct((plink*) new_v);

#ifndef NoShortcut
      // only repeat if oldv was indirect
      // if(ii == 0 && (oldv == nullptr || !oldv->is_indirect())) break;

    }
#endif
    V* curv = v.load();
    if(curv != nullptr) set_stamp(curv);
    return false;
  }

  V* operator=(V* b) {store(b); return b; }
};
} // namespace vl
