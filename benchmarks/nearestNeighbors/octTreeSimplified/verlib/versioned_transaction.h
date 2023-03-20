#pragma once

#include "flock/flock.h"
#include "flock/acquired_pool.h"
#include "timestamps.h"
#include "transaction_descriptor.h"
#include "memory_pool.h"

namespace verlib {

  bool valid_pointer(void* ptr) { return ((size_t) ptr) < (1ul << 48);}

  // just used for correctness assertions
#define init_ptr (((versioned*) nullptr) - 1)

  struct versioned {
  private:
    flck::atomic_write_once<TS> stamp_or_descriptor;
    flck::atomic_write_once<versioned*> next_version;

    // tags second highest bit to mark as descriptor
    constexpr static int bit_num = 62;
    static bool is_tagged(TS sd) {return (sd >> bit_num) & 1;}
    static TS tag_descriptor(trans_descriptor* ptr) {
      return (TS) ptr | (1ul << bit_num);}
    static trans_descriptor* untag_descriptor(TS sd) {
      auto ptr = (trans_descriptor*) (sd & ~(1ul << bit_num));
      assert(valid_pointer((void*) ptr));
      return ptr;
    }

  public:

    std::pair<trans_descriptor*, TS> get_descriptor_and_stamp() {
      TS td = stamp_or_descriptor.load();
      if (is_tagged(td)) {
      	trans_descriptor* descriptor = untag_descriptor(td);
      	trans_descriptor_pool.acquire(descriptor);
      	td = stamp_or_descriptor.load();
      } 
      if (is_tagged(td)) return std::pair(untag_descriptor(td), tbd);
      else return std::pair(nullptr, td);
    }

    void set_descriptor(trans_descriptor* descriptor) {
      assert(stamp_or_descriptor.load() == tbd);
      stamp_or_descriptor.init(tag_descriptor(descriptor));
    }

    void remove_descriptor(TS ts, bool check_stamp = true) {
#ifndef NDEBUG
      TS td = stamp_or_descriptor.load();
#endif
      assert(is_tagged(td));
      assert(untag_descriptor(td)->owner == flck::internal::current_id);
      assert(!check_stamp || untag_descriptor(td)->time_stamp.load() == ts);
      stamp_or_descriptor = ts;
    }

    // once a next version is added, it becomes recorded
    bool is_recorded() { return next_version.load() != init_ptr; }

    void set_stamp(TS t) {
      assert(stamp_or_descriptor.load() == tbd);
      next_version = nullptr;
      stamp_or_descriptor = t;
    }

    // use when stamp is already set
    TS get_stamp() {
      TS td = stamp_or_descriptor.load();
      assert(!is_tagged(td));
      return td;
    }

    // use when stamp could be descriptor
    TS get_stamp_indirect() {
      TS td = stamp_or_descriptor.load();
      if (is_tagged(td)) return tbd;
      else return td;
    }

    versioned() : stamp_or_descriptor(tbd), next_version(init_ptr) {} 
    versioned(versioned* next) : stamp_or_descriptor(tbd), next_version(next) {}

    versioned* get_next_version() {return next_version.load();}
    void set_next_version(versioned* nv) { next_version = nv;}
    void init_next_version(versioned* nv) {
      assert(next_version.load() == init_ptr);  
      next_version.init(nv);}

    void set_with_current_stamp() {
      TS old_t = tbd;
      TS new_t = global_stamp.get_write_stamp();
      TS td = stamp_or_descriptor.load_ni();
      assert(!is_tagged(td));    
      if (td == old_t) 
	stamp_or_descriptor.cas_ni(old_t, new_t);
    }
  };

  // uses lowest bit of pointer to indicate whether indirect (1) or not (0)
  versioned* add_indirect(versioned* ptr) {
    return (versioned*) (1ul | (size_t) ptr);};
  versioned* strip_indirect(versioned* ptr) {
    return (versioned*) ((size_t) ptr & ~1ul);}
  bool is_indirect(versioned* ptr) {
    return (size_t) ptr & 1;}

  struct ver_link : versioned {
    versioned* value;
    ver_link(versioned* next, versioned* value) : versioned{next}, value(value) {}
  };

  verlib::memory_pool<ver_link> link_pool;

  void shortcut(flck::atomic<versioned*>* loc, versioned* ptr) {
    ver_link* ptr_ = (ver_link*) strip_indirect(ptr);
    flck::non_idempotent([=] {
      if (ptr_->get_stamp_indirect() <= done_stamp) 
	if (loc->cas(ptr, ptr_->value)) link_pool.pool.retire(ptr_);});
  }

  // Installs ptr to front of version list for a store on loc.
  // - If in a transaction, then set descriptor on the new link, and the
  // store will be completed by finalize_versioned_store.
  // - Otherwise set the timestamp and no finalization is required.
  void install_versioned_store(flck::atomic<versioned*>* loc,
			       versioned* ptr,
			       trans_descriptor* descriptor) {
    assert(current_transaction == nullptr);
    versioned* old_v = loc->load();
    versioned* new_v = ptr;
    // indirection on need
    bool use_indirect = (ptr == nullptr || ptr->is_recorded()); //log 
    if (use_indirect) {
      ver_link* link = link_pool.new_obj(old_v, new_v);
      link_pool.retire_on_abort(descriptor, link);
      new_v = add_indirect(link);
    }
    else ptr->init_next_version(old_v);

    // install descriptor if in transaction
    if (descriptor != nullptr)
      strip_indirect(new_v)->set_descriptor(descriptor); // could use an init
#ifdef NoShortcut
    *loc = new_v;
    if (is_indirect(old_v))
      link_pool.retire_on_success(descriptor, (ver_link*) strip_indirect(old_v)); 
#else
    // if using shortcutting then need to deal with race
    loc->cam(old_v, new_v);  // load is logged
    if (is_indirect(old_v)) {
      versioned* val = loc->load(); // logged
      if (val == new_v)
	link_pool.retire_on_success(descriptor, (ver_link*) strip_indirect(old_v)); 
      else {
	strip_indirect(new_v)->set_next_version(val); // logged, but unusual
	loc->cam(val, new_v); // logged
      }
    }
#endif
    // set timestamp and try shortcut if not in transaction
    if (descriptor == nullptr) {
      strip_indirect(new_v)->set_with_current_stamp(); // need not be idempotent since self set
      if (use_indirect) shortcut(loc, new_v); // need not be idempotent
    }
    //return old_v;
  }

  void finalize_versioned_store(flck::atomic<versioned*>* loc,
				TS time_stamp) {
    versioned* ptr = strip_indirect(loc->load()); // could avoid load by storing pointer when installed
    assert(ptr != nullptr);
    // splice this latest version out if transaction aborted
    if (time_stamp == failed) *loc = ptr->get_next_version();  // need not be logged
    // replace descriptor with timestamp (possibly failed stamp)
    ptr->remove_descriptor(time_stamp); // write once, so no commit
  }

  void validate_reads(trans_descriptor* descriptor);

  // If a descriptor is installed then need to resolve whether to
  // return the new value (if succeded) or old value (if not yet in
  // validation stage or if failed).
  template <typename T>
  T resolve_descriptor(trans_descriptor* descriptor, T old_v, T new_v,
		       bool validating) {
    if (descriptor->owner == flck::internal::current_id) {
      assert(validating); // self descriptor should not be installed if not validating
      return old_v; // if validating then should check against old value
    }
    TS desc_ts = descriptor->time_stamp.load();
    if (desc_ts == tbd || desc_ts == failed) return old_v;
#ifdef Validate
    if (validating && desc_ts == tbs &&
	descriptor->owner < flck::internal::current_id) { 
      // kill other if they have lower id. prevents validate cycle
      TS old = tbs;
      descriptor->time_stamp.cam(old, failed);
    } else { //otherwise wait for a bit to see if other is done, then help
      int cnt = 0;
      while (descriptor->time_stamp.load() == tbs && cnt++ < 50)
      	for (volatile int i=0; i < 100; i++);
      if (descriptor->time_stamp.load() == tbs) validate_reads(descriptor);
    }
    if (descriptor->time_stamp.load() == failed) return old_v;
#else
    descriptor->atomic_set_current_stamp();
#endif
    return new_v;
  }

  versioned* get_value(versioned* tagged_ptr, bool validating = false) {
    if (tagged_ptr == nullptr) return nullptr;
    versioned* ptr = strip_indirect(tagged_ptr);
    auto [descriptor, stamp] = ptr->get_descriptor_and_stamp();
    if (stamp == failed) return ptr->get_next_version(); // failed so ignore value
    if (stamp != tbd) return tagged_ptr; // value already fully set
    if (descriptor == nullptr) { // store was not in a transaction
      ptr->set_with_current_stamp(); // check and set the stamp
      return tagged_ptr; }
    // a transaction descriptor is installed, need to resolve the state
    return resolve_descriptor(descriptor, ptr->get_next_version(),
			      tagged_ptr, validating);
  }

  bool validate_pointer(flck::atomic<versioned*>* location, versioned* expected) {
    versioned* ptr = get_value(location->load(), true);
    // second condition is to deal with shortcutting
    return (ptr == expected || (is_indirect(expected) &&
				ptr == ((ver_link*) strip_indirect(expected))->value));
  }

  thread_local int try_count;
  thread_local bool trans_aborted = false;
#include "lock.h"
#include "transactions.h"
  //  #include "atomic_val.h"

  template <typename V>
  struct versioned_ptr {
  private:
    flck::atomic<versioned*> v;

    // static V* get_ptr(versioned* ptr) {
    //   versioned* ptr_ = strip_indirect(ptr);
    //   if (!is_indirect(ptr)) return (V*) ptr_;
    //   return (V*) ((ver_link*) ptr_)->value;
    // }

    V* get_ptr_shortcut(versioned* ptr) {
      versioned* ptr_ = strip_indirect(ptr);
      if (!is_indirect(ptr)) return (V*) ptr_;
#ifndef NoShortcut
      shortcut(&v, ptr);
#endif
      return (V*) ((ver_link*) ptr_)->value;
    }

    static versioned* set_zero_stamp(V* ptr) {
      if (ptr != nullptr && !ptr->is_recorded()) 
	ptr->set_stamp(zero_stamp);
      return ptr;
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
      versioned* ptr = get_value(v.read());
      versioned* ptr_ = strip_indirect(ptr);

      // chase down version chain
      while (ptr != nullptr && ptr_->get_stamp() > ls) {
	ptr = ptr_->get_next_version();
	ptr_ = strip_indirect(ptr);
      }
#ifdef LazyStamp
      if (ptr != nullptr && ptr_->get_stamp() == ls
	  && speculative)
	aborted = true;
#endif
      if (is_indirect(ptr)) {
	return (V*) ((ver_link*) ptr_)->value;
      } else return (V*) ptr;
    }

    std::optional<int> find_in_write_log() {
      auto& log = current_transaction->write_ptr_log;
      for (int i=0; i < log.size(); i++)
	if (log[i].first == &v)
	  return std::optional<int>(i);
      return std::optional<int>();
    }
      
    V* load() {  
      assert(this != nullptr);
      if (local_stamp != -1) return read_snapshot();
      if (current_transaction != nullptr) {
	// first see if it is in the write log
	std::optional<int> lid = find_in_write_log();
	if (lid.has_value())
	  return (V*) current_transaction->write_ptr_log[*lid].second;

	// if not then get the value and add to validate log
	versioned* cur = get_value(v.load());
	//if (in_lock) validate_ptr_log.add(std::tuple(&v, cur));
	if (in_lock) current_transaction->validate_ptr_log.add(std::tuple(&v, cur));

	// for opacity check that has not changed since start of transaction
	if (cur != nullptr && strip_indirect(cur)->get_stamp_indirect() > start_stamp)
	  trans_aborted = true;
	return get_ptr_shortcut(cur);
      }
      return get_ptr_shortcut(get_value(v.load()));
    }

    void store(V* ptr) { 
      if (current_transaction != nullptr) {
	// check if it is in the write log: update if so, otherwise add
	std::optional<int> lid = find_in_write_log();
	if (lid.has_value()) current_transaction->write_ptr_log[*lid].second = ptr;
	else current_transaction->write_ptr_log.add(std::pair(&v,ptr));
      } else { // not in a transaction
	install_versioned_store(&v, ptr, nullptr);
      }
    }
    V* operator=(V* b) {store(b); return b; }

    bool validate(V* expected) {
      assert(this != nullptr);
      // check if in write log.  if so no need to validate
      std::optional<int> lid = find_in_write_log();
      if (lid.has_value())
	return (current_transaction->write_ptr_log[*lid].second == expected);

      // otherwise check that value matches expected, and if so add to validate log
      versioned* ptr = get_value(v.load());
      V* cur = (V*) (is_indirect(ptr) ? ((ver_link*) strip_indirect(ptr))->value : ptr);
      if (cur != expected) return false;
      // Note that adding "ptr" not "expected" to log. This is because
      // ptr might be an indirect pointer to the expected value and
      // the pointers need to match to avoid ABA problems.
      //validate_ptr_log.add(std::tuple(&v, ptr));
      current_transaction->validate_ptr_log.add(std::tuple(&v, ptr));
      return true;
    }
  };

  struct Empty : versioned {};
  Empty true_val;

  // Simulates a boolean with a pointer to an "Empty" dummy struct.
  // If false, then pointer is nullptr otherwise poiter to "true_val".
  // Avoids level of indirection when storing false.
  struct atomic_bool {
    versioned_ptr<Empty> val;
    bool load() {return val.load() != nullptr;}
    void store(bool v) {val = v ? &true_val : nullptr;}
    bool validate(bool expected) {
      Empty* cur = val.load();
      if ((cur == nullptr) == expected) return false;
      return val.validate(cur); }
    atomic_bool(bool v) : val(v ? &true_val : nullptr) {}
    bool operator=(bool b) {store(b); return b; }
  };

  // should probably go into the lock code.
  template <typename T, typename E>
  bool validate(verlib::lock& lck, T* v, E expected) {
    if (current_transaction != nullptr) {
#ifdef Validate
      if (true) //try_count < 100) // anywhere from 10-100 seems good
	return v->validate(expected);
      else // if too many tries give up on validation and take locks
	return lck.try_lock([=] {return (v->load() == expected);});
#else
      return lck.try_lock([=] {return (v->load() == expected);});
#endif
    } else return true;
  }


} // namespace verlib
