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

namespace verlib {
  
#define bad_ptr ((void*) ((1ul << 48) -1))

struct versioned {
  std::atomic<TS> time_stamp;
  void* next_version;
  versioned() : time_stamp(tbd), next_version(bad_ptr) {}
};

template <typename V>
struct versioned_ptr {
private:
  std::atomic<V*> v;
  // sets the timestamp in a version link if time stamp is TBD
  static V* set_stamp(V* x) {
    assert(x != nullptr);
    if (x->time_stamp.load() == tbd) {
      TS ts = global_stamp.get_write_stamp();
      long old = tbd;
      if (x->time_stamp.load() == tbd)
        x->time_stamp.compare_exchange_strong(old, ts);
    }
    return x;
  }

  static V* set_zero(V* x) {
    if (x != nullptr && x->time_stamp.load() == tbd)
      x->time_stamp = zero_stamp;
    return x;
  }

public:

  versioned_ptr(V* v) : v(set_zero(v)) {}
  versioned_ptr(): v(nullptr) {}
  void init(V* vv) { v = set_zero(vv);}

  // reads snapshotted version (ls >= 0)
  V* read_snapshot() {
    // ensure time stamp is set
    V* head = v.load();
    if(head == nullptr) return nullptr;
    set_stamp(head);
    TS ls = local_stamp;
    while (head->time_stamp.load() > ls) {
      head = (V*) head->next_version;
    }
#ifdef LazyStamp
    if (head->time_stamp.load() == ls && speculative)
      aborted = true;
#endif
    return head;
  }

  V* load() {
    if (local_stamp != -1) return read_snapshot();
    else {
      V* head = flck::commit(v.load());
      if(head != nullptr) set_stamp(head);
      return head;
    }
  }

  V* read() {
    if (local_stamp != -1) return read_snapshot();
    return v.load();}

  V* read_cur() {
    return v.load();}

  void validate() { 
    V* head = v.load();
    if(head != nullptr) set_stamp(head);
  }

#ifdef SlowStore
  void store(V* newv) {
    V* oldv = load();
    cas(oldv, newv);
  }
#else  
  void store(V* newv) {
    V* oldv = flck::commit(v.load());
    if (newv == nullptr) {
      std::cout << "recording with nullptr not allowed" << std::endl;
      abort();
    }
    // check that newv is only recorded once
    if (flck::commit(newv->time_stamp.load()) != tbd) {
      std::cout << "recording a second time not allowed" << std::endl;
      abort();
    }

    flck::skip_if_done_no_log([&] { // for efficiency, correct without it
      newv->next_version = (void*) oldv;
      v.compare_exchange_strong(oldv, newv);
      set_stamp(newv);
      // shortcut if appropriate
      //if (oldv != nullptr && newv->time_stamp == oldv->time_stamp)
      //newv->next_version = oldv->next_version;
    });
  }
#endif

  bool cas(V* expv, V* newv) {
    if (newv == nullptr) {
      std::cout << "recording with nullptr not allowed" << std::endl;
      abort();
    }

    V* oldv = v.load();
    if(oldv != nullptr) set_stamp(oldv);
    if(oldv != expv) return false;
    if(oldv == newv) return true;
    newv->next_version = expv;
    if(v.compare_exchange_strong(oldv, newv)) {
      set_stamp(newv);
      return false;
    } else {
      set_stamp(v.load());
      return true;
    }
  }

  V* operator=(V* b) {store(b); return b; }
};

} // namespace verlib
