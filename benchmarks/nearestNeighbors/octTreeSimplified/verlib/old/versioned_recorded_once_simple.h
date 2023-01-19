// A version for data structures that only record once.
// Effectively this means that an object that is pointed to
// by a versioned_ptr can only be stored to a versioned pointer
// once.

// Based on paper:
// Yuanhao Wei, Naama Ben-David, Guy E. Blelloch, Panagiota Fatourou, Eric Ruppert, Yihan Sun
// Constant-time snapshots with applications to concurrent data structures
// PPoPP 2021
#pragma once
#include "flock/flock.h"
#include "timestamps.h"

struct versioned {
  std::atomic<TS> time_stamp;
  versioned* next_version;
  versioned() : time_stamp(tbd) {}
};

template <typename V>
struct versioned_ptr {
private:
  flck::atomic<V*> v;

  static V* set_stamp(V* ptr) {
    if (ptr->time_stamp.load() == tbd) {
      TS old_t = tbd;
      TS new_t = global_stamp.get_write_stamp();
      ptr->time_stamp.compare_exchange_strong(old_t, new_t);
    }
    return ptr;
  }

  static V* init_ptr(V* x) {
    if (x->time_stamp.load() == tbd)
      x->time_stamp = zero_stamp;
    return x;
  }

public:

  versioned_ptr(V* v) : v(init_ptr(v)) {}
  versioned_ptr(): v(nullptr) {}
  void init(V* vv) { v = init_ptr(vv);}

  V* read_snapshot() {
    V* head = set_stamp(v.load());
    while (head->time_stamp.load() > local_stamp)
      head = (V*) head->next_version;
#ifdef LazyStamp
    if (head->time_stamp.load() == local_stamp)
      bad_stamp = true;
#endif
    return head;
  }

  V* load() {
    if (local_stamp != -1) return read_snapshot();
    else return set_stamp(v.load());}

  V* read() { return v.read();}

  void validate() { set_stamp(v.load());}
  
  void store(V* new_v) {
    V* old_v = v.load();
    new_v->next_version = (versioned*) old_v;
    v = new_v;
    set_stamp(new_v);
  }

  bool cas(V* exp_v, V* new_v) {
    V* old_v = v.load();
    if (old_v != nullptr) set_stamp(old_v);
    if (old_v != exp_v) return false;
    if (old_v == new_v) return true;
    new_v->next_version = exp_v;
    if (v.cas_ni(old_v, new_v)) {
      set_stamp(new_v);
      return true;
    } else {
      set_stamp(v.load());
      return true;
    }
  }

  V* operator=(V* b) {store(b); return b; }
};
