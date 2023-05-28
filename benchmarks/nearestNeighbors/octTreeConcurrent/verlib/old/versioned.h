#pragma once
#include "flock/flock.h"
#include "timestamps.h"

namespace vl {
#define bad_ptr ((1ul << 48) -1)

using IT = size_t;
struct versioned {
  std::atomic<TS> time_stamp;
  std::atomic<IT> next_version;
  versioned() : time_stamp(tbd), next_version(bad_ptr) {}
  versioned(IT next) : time_stamp(tbd), next_version(next) {}
};

struct plink : versioned {
  IT value;
  plink(IT next, IT value) : versioned(next), value(value) {}  
};

 flck::internal::mem_pool<plink> link_pool;

template <typename V>
struct versioned_ptr {
private:
  using TV = flck::internal::tagged<V*>;
  std::atomic<IT> v;

  // uses lowest three bits as mark:
  //   2nd bit to indicate it is an indirect pointer
  //   1st bit to indicate it is a null pointer via an indirect pointer (2nd bit also set)
  //   3rd bit to indicate time_stamp has not been set yet
  // the highest 16 bits are used by the ABA "tag"
  static V* add_null_mark(V* ptr) {return (V*) (3ul | (IT) ptr);};
  static V* add_indirect_mark(V* ptr) {return (V*) (2ul | (IT) ptr);};
  static bool is_empty(IT ptr) {return ptr & 1;}
  static bool is_indirect(IT ptr) {return (ptr >> 1) & 1;}
  static bool is_null(IT ptr) {return TV::value(ptr) == nullptr || is_empty(ptr);}
  static V* strip_mark_and_tag(IT ptr) {return TV::value(ptr & ~7ul);}
  static V* get_ptr(IT ptr) {
    return (is_indirect(ptr) ?
	    (is_empty(ptr) ? nullptr : (V*) ((plink*) strip_mark_and_tag(ptr))->value) :
	    strip_mark_and_tag(ptr));}
  

  // sets the timestamp in a version link if time stamp is TBD
  // The test for is_null is an optimization avoiding reading the timestamp
  // for indirectly stored nullptrs
  static IT set_stamp(IT newv) {
    if (!is_null(newv)) {
      V* x = strip_mark_and_tag(newv);
      if ((x != nullptr) && x->time_stamp.load() == tbd) {
      	TS ts = global_stamp.get_write_stamp();
      	long old = tbd;
        if(x->time_stamp.load() == tbd)
      	  x->time_stamp.compare_exchange_strong(old, ts);
      }
    }
    return newv;
  }

  static V* set_zero(V* ptr) {
    if (ptr != nullptr && ptr->time_stamp.load() == tbd)
      ptr->time_stamp = zero_stamp;
    return ptr;
  }

  // For an indirect pointer if its stamp is older than done_stamp
  // then it will no longer be accessed and can be spliced out.
  std::pair<V*,bool> shortcut_indirect(IT ptr) {
    auto ptr_notag = (plink*) strip_mark_and_tag(ptr);
    if (is_indirect(ptr)) {
      //i_counts[16*parlay::worker_id()]++;
      TS stamp = ptr_notag->time_stamp.load();
      V* newv = (V*) ptr_notag->value;
      if (stamp <= done_stamp) {
        // there can't be an ABA unless indirect node is reclaimed
        if (TV::cas_with_same_tag(v, ptr, newv, true)) 
          link_pool.retire(ptr_notag);
        return std::pair{newv, true};
      } else return std::pair{newv, false};
    } else return std::pair{(V*) ptr_notag, false};
  }

public:

  versioned_ptr(V* v) : v(TV::init(set_zero(v))) {}
  versioned_ptr(): v(TV::init(0)) {}
  ~versioned_ptr() {
    IT ptr = v.load();
    if (is_indirect(ptr))
      link_pool.destruct((plink*) strip_mark_and_tag(ptr));
  }

  void init(V* vv) {v = TV::init(set_zero(vv));}
  // reads snapshotted version
  V* read_snapshot() {
    TS ls = local_stamp;
    IT head = v.load();
    set_stamp(head);
    V* head_unmarked = strip_mark_and_tag(head);

    // chase down version chain
    while (head_unmarked != 0 && head_unmarked->time_stamp.load() > ls) {
      head = head_unmarked->next_version.load();
      head_unmarked = strip_mark_and_tag(head);
    }
#ifdef LazyStamp
    if (head_unmarked != 0 && head_unmarked->time_stamp.load() == ls && speculative)
      aborted = true;
#endif
    return get_ptr(head);
  }

  V* load() {  // can be used anywhere
    if (local_stamp != -1) return read_snapshot();
    else return shortcut_indirect(set_stamp(flck::commit(v.load()))).first;}

  V* read() {  // only safe on journey
    //if (local_stamp != -1) read_snapshot();
    return shortcut_indirect(v.load()).first;
  }

  V* read_cur() {  // only safe on journey, outside of snapshot
    return shortcut_indirect(v.load()).first;
  }

  void validate() {
    set_stamp(v.load());     // ensure time stamp is set
  }
  
  void store(V* newv_) {
    IT oldv_tagged = flck::commit(v.load());
    V* newv = newv_;
    V* newv_marked = newv;

    flck::skip_if_done([&] {
      if (newv_ == nullptr || flck::commit(newv_->time_stamp.load() != tbd)) {
        newv = (V*) link_pool.new_obj((IT) oldv_tagged, (IT) newv);
        newv_marked = add_indirect_mark(newv);
      } else {
        IT initial_ptr = bad_ptr;
        if(newv->next_version == initial_ptr)
          newv->next_version.compare_exchange_strong(initial_ptr, 
            (IT) oldv_tagged);
      }
      bool succeeded = TV::cas(v, oldv_tagged, newv_marked);
      IT x = flck::commit(v.load());
      
      if (is_indirect(oldv_tagged)) {
        if (succeeded) link_pool.retire((plink*) strip_mark_and_tag(oldv_tagged));
        else if (TV::get_tag(x) == TV::get_tag(oldv_tagged)) { 
          succeeded = TV::cas(v, x, newv_marked);
         	x = flck::commit(v.load());
        }
      } 

      // now set the stamp from tbd to a real stamp
      set_stamp(x);

      // try to shortcut indirection out, and if not, clear unset mark
      // for time stamp
      if (!shortcut_indirect(x).second && newv_ == nullptr)
        TV::cas(v, x, add_null_mark(TV::value(x)));  // clear the "unset" mark

      // shortcut version list if appropriate, getting rid of redundant
      // time stamps.  
      // if (oldv != nullptr && newv->time_stamp == oldv->time_stamp) 
      //   newv->next_version = oldv->next_version.load();

      // free if allocated link was not used
      if (!succeeded && is_indirect((IT) newv_marked))
        link_pool.destruct((plink*) newv);
    });
  }
  
  bool cas(V* expv_, V* newv_) {
    // CAS could be interuupted once by shortcut indirect and once by setting null bit
    for(int ii = 0; ii < 3; ii++) {
      IT oldv_tagged = v.load();
      set_stamp(oldv_tagged);
      if(shortcut_indirect(oldv_tagged).first != expv_) return false;
      if(expv_ == newv_) return true;

      V* newv = newv_;
      V* newv_marked = newv;

      bool use_indirect = (newv_ == nullptr || newv_->time_stamp.load() != tbd);

      if(use_indirect) {
        newv = (V*) link_pool.new_obj((IT) oldv_tagged, (IT) newv);
        newv_marked = add_indirect_mark(newv);
      } else newv->next_version = (IT) oldv_tagged;

      // swap in new pointer
      IT x = oldv_tagged;
      // use cas_with_same_tag to avoid picking new timestamp
      bool succeeded = TV::cas_with_same_tag(v, x, newv_marked); 
      
      if(succeeded) {
        set_stamp((IT) newv_marked);
        // try to shortcut indirection out, and if not, clear unset mark
        // for time stamp
        if(is_indirect(oldv_tagged))
          link_pool.retire((plink*) strip_mark_and_tag(oldv_tagged));
        if (use_indirect && !shortcut_indirect((IT) newv_marked).second && newv_ == nullptr)
          TV::cas_with_same_tag(v, (IT) newv_marked, add_null_mark(newv_marked)); 
        return true;
      }
      if(use_indirect) link_pool.destruct((plink*) newv);
    }

    set_stamp((IT) strip_mark_and_tag(v.load()));
    return false;

    // // shortcut version list if appropriate, getting rid of redundant
    // // time stamps.  
    // if (oldv != nullptr && newv->time_stamp == oldv->time_stamp) 
    //   newv->next_version = oldv->next_version.load();
  }

  V* operator=(V* b) {store(b); return b; }
};
} // namespace verlib
