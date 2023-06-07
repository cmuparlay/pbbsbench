#pragma once

#include "flock/flock.h"
#include "timestamps.h"

namespace verlib {
  using flck::memory_pool;

const TS tbd = std::numeric_limits<TS>::max()/4;

  template <typename F>
  auto do_now(F f) {return f();}

template <typename T>
using atomic = flck::atomic<T>;
using lock = flck::lock;
using atomic_bool = flck::atomic_write_once<bool>;
  
struct versioned {
  flck::atomic_write_once<TS> time_stamp;
  versioned* next_version;
  versioned() : time_stamp(tbd) {}
  versioned(versioned* next) : time_stamp(tbd), next_version(next) {}
};

struct ver_link : versioned {
  versioned* value;
  ver_link(versioned* next, versioned* value) : versioned{next}, value(value) {}  
};

flck::memory_pool<ver_link> link_pool;

template <typename V>
struct versioned_ptr {
private:
  flck::atomic<versioned*> v;

  // uses lowest bit of pointer to indicate whether indirect (1) or not (0)
  static versioned* add_indirect(versioned* ptr) {
    return (versioned*) (1ul | (size_t) ptr);};
  static versioned* strip_indirect(versioned* ptr) {
    return (versioned*) ((size_t) ptr & ~1ul);}
  static bool is_indirect(versioned* ptr) {
    return (size_t) ptr & 1;}

  void shortcut(versioned* ptr) {
#ifndef NoShortcut
    ver_link* ptr_ = (ver_link*) strip_indirect(ptr);
    if (ptr_->time_stamp.load_ni() <= done_stamp) {
#ifdef NoHelp
      if (v.cas(ptr, ptr_->value))
	link_pool.retire(ptr_);
#else
      if (v.cas_ni(ptr, ptr_->value))
	link_pool.retire_ni(ptr_);
#endif
    }
#endif
  }

  static V* get_ptr(versioned* ptr) {
    versioned* ptr_ = strip_indirect(ptr);
    if (!is_indirect(ptr)) return (V*) ptr_;
    return (V*) ((ver_link*) ptr_)->value;
  }

  V* get_ptr_shortcut(versioned* ptr) {
    versioned* ptr_ = strip_indirect(ptr);
    if (!is_indirect(ptr)) return (V*) ptr_;
    shortcut(ptr);
    return (V*) ((ver_link*) ptr_)->value;
  }

  static versioned* set_stamp(versioned* ptr) {
    versioned* ptr_ = strip_indirect(ptr);
    if (ptr != nullptr && ptr_->time_stamp.load_ni() == tbd) {
      TS t = global_stamp.get_write_stamp();
      if(ptr_->time_stamp.load_ni() == tbd)
        ptr_->time_stamp.cas_ni(tbd, t);
    }
    return ptr;
  }

  static versioned* set_zero_stamp(V* ptr) {
    if (ptr != nullptr && ptr->time_stamp.load_ni() == tbd)
      ptr->time_stamp = zero_stamp;
    return ptr;
  }

  bool cas_from_cam(versioned* old_v, versioned* new_v) {
#ifdef NoHelp
    return v.cas(old_v, new_v);
#else
    v.cam(old_v, new_v);
    return (v.load() == new_v ||
	    strip_indirect(new_v)->time_stamp.load() != tbd);
#endif
  }

public:

  versioned_ptr(): v(nullptr) {}
  versioned_ptr(V* ptr) : v(set_zero_stamp(ptr)) {}

  ~versioned_ptr() {
    versioned* ptr = v.load();
    if (is_indirect(ptr))
      link_pool.destruct((ver_link*) strip_indirect(ptr));
  }

  void init(V* ptr) {v = set_zero_stamp(ptr);}

  V* read_snapshot() {
    TS ls = local_stamp;
    versioned* head = set_stamp(v.read());
    versioned* head_unmarked = strip_indirect(head);

    // chase down version chain
    while (head != nullptr && head_unmarked->time_stamp.load() > ls) {
      head = head_unmarked->next_version;
      head_unmarked = strip_indirect(head);
    }
#ifdef LazyStamp
    if (head != nullptr && head_unmarked->time_stamp.load() == ls
	&& speculative)
      aborted = true;
#endif
    if (is_indirect(head)) {
      return (V*) ((ver_link*) head_unmarked)->value;
    } else return (V*) head;
  }

  V* load() {  // can be used anywhere
    if (local_stamp != -1) return read_snapshot();
    else return get_ptr_shortcut(set_stamp(v.load()));
  }
  
  V* read() {  // only safe on journey
    return get_ptr_shortcut(v.read());
  }

  void validate() {
    set_stamp(v.load());     // ensure time stamp is set
  }

#ifdef SlowStore
  void store(V* newv) {
    V* oldv = load();
    cas(oldv, newv);
  }
#else  
  void store(V* ptr) {
    versioned* old_v = v.load();
    versioned* new_v = ptr;
    bool use_indirect = (ptr == nullptr || ptr->time_stamp.load() != tbd);

    if (use_indirect) 
      new_v = add_indirect(link_pool.new_obj(old_v, new_v));
    else ptr->next_version = old_v;

#ifdef NoShortcut
    v = new_v;
    if (is_indirect(old_v))
      link_pool.retire((ver_link*) strip_indirect(old_v));
#else
    v.cam(old_v, new_v);
    if (is_indirect(old_v)) {
      versioned* val = v.load();
      ver_link* old_l = (ver_link*) strip_indirect(old_v);
      if (val != old_l->value)
	link_pool.retire(old_l);
      else v.cam(val, new_v);
    }
#endif
    set_stamp(new_v);
    if (use_indirect) shortcut(new_v);
  }
#endif
  
  bool cas(V* exp, V* ptr) {
#ifndef NoShortcut
    for(int ii = 0; ii < 2; ii++) {
#endif
      versioned* old_v = v.load();
      versioned* new_v = ptr;
      V* old = get_ptr(old_v);
      set_stamp(old_v);
      if(old != exp) return false;
      if (exp == ptr) return true;
      bool use_indirect = (ptr == nullptr || ptr->time_stamp.load() != tbd);

      if(use_indirect)
	new_v = add_indirect(link_pool.new_obj(old_v, new_v));
      else ptr->next_version = old_v;

      if (cas_from_cam(old_v, new_v)) {
        set_stamp(new_v);
        if (is_indirect(old_v))
	  link_pool.retire((ver_link*) strip_indirect(old_v));
#ifndef NoShortcut
	if (use_indirect) shortcut(new_v);
#endif
        return true;
      }
      if (use_indirect)
	link_pool.destruct((ver_link*) strip_indirect(new_v));
#ifndef NoShortcut
    }
#endif
    set_stamp(v.load());
    return false;
  }

  V* operator=(V* b) {store(b); return b; }
};

  template <typename T, typename E>
  bool validate(flck::lock& lck, T* v, E expected) {
    return true;
  }

} // namespace verlib
