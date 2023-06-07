#pragma once
#include <atomic>
#include "epoch.h"

namespace flck {
namespace internal {

// a flag to indicate that currently helping
static thread_local bool helping = false;

// default log length.  Will grow if needed.
constexpr int Log_Len = 8;
using log_entry = std::atomic<void*>;
struct log_array;

mem_pool<log_array> log_array_pool;

struct log_array {
  std::array<log_entry,Log_Len> log_entries;
  std::atomic<log_array*> next;

  void init() { // TODO: can this just be a constructor?
    for (int i=0; i < Log_Len; i++) {
      if (log_entries[i].load() != nullptr) 
        log_entries[i].store(nullptr, std::memory_order::memory_order_relaxed);
    }
    if(next.load() != nullptr) 
      next.store(nullptr, std::memory_order::memory_order_relaxed);
  }

  log_entry & operator [](int i) { return log_entries[i]; }

  ~log_array() {
    if(next != nullptr) {
      log_array_pool.destruct(next);
    }
  }
};

// A log_array is associated with a thunk and keeps track of
// "important events" that all threads helping the thunk need to agree
// on.  Each thread maintains a log, which is a pointer into the
// common log_array, and the current position (count) it is in the
// log_array (initially 0).  Each thread tries to commit the important
// events to the shared log_array and whoever does so first succeeds.
// The others read the results of the commited value to agree on the
// value.  After each commit a thread advances its pointer to the next
// location in the log array.  The log_array consists of a fixed size
// array of the log entries and a next pointer (initially null).  If a
// log array runs out of slots a new log is allocated and is linked to
// by the next pointer.  Hence logs have unbounded length.  A log
// entry of zero means it has not been committed.
//
// commit_value commits to the log and advances the count.  Its value
// cannot be zero.  More details below.
//
// commit_value_safe allows the value to be zero but uses bit 48 as a
// marked flag so it only supports up to 6 bytes It supports pointers
// assuming they fit within 6 bytes (i.e. upper 2 bytes are empty).
struct Log {
  log_array* vals;
  int count;
  Log(log_array* pa, int c) : vals(pa), count(c) {}
  Log() : vals(nullptr), count(0) {}

  // get the next entry in the log
  log_entry* next_entry() {
    assert(!is_empty());
    if (count == Log_Len) { // ran out of slots
      count = 0;
      log_array* next_log_array = vals->next;
      if (next_log_array != nullptr) vals = next_log_array;
      else {  // next_log_array == nullptr, try to commit a new log array
        log_array* new_log_array = log_array_pool.new_obj();
        new_log_array->init();
        if(vals->next.compare_exchange_strong(next_log_array, new_log_array))
          vals = new_log_array;
        else { // if in meantime someone else allocated it, then return memory
          vals = next_log_array;
          log_array_pool.destruct(new_log_array);
        }
      }
    }
    return &(*vals)[count++];
  }
  log_entry* current_entry() {return &(*vals)[count-1];}
  bool is_empty() {return vals == nullptr;}

  // commits a value to the log, or returns existing value if already committed
  // along with a false flag.
  // V must be convertible to void*
  // initialized with null pointer, so should not commit a nullptr (or zero)
  // code tags pointers with a count or flag to avoid this
  template<typename V>
  std::pair<V,bool> commit_value(V newv) {
    if (is_empty()) return std::make_pair(newv, true);
    assert(newv != (V) 0); // check not committing 0
    log_entry* l = next_entry();
    void* oldv = l->load();
    if (oldv == nullptr && l->compare_exchange_strong(oldv, (void*) newv))
      return std::make_pair(newv, true);
    else return std::make_pair((V) (size_t) oldv, false);
  }

  // this version tags 48th bit so value can be zero
  template<typename V>
  std::pair<V,bool> commit_value_safe(V val) {
    static_assert(sizeof(V) <= 6 || std::is_pointer<V>::value ||
		  "Type for commit_value_safe must be a pointer or at most 6 bytes");
    if (is_empty()) return std::make_pair(val, true);
    size_t set_bit = (1ul << 48);
    log_entry* l = next_entry();
    void* oldv = l->load();
    void* newv = (void*) (((size_t) val) | set_bit);
    if (oldv == nullptr && l->compare_exchange_strong(oldv, (void*) newv))
      return std::make_pair(val, true);
    //else return std::make_pair((V) (((size_t) oldv) & (set_bit - 1)), false);
    else return std::make_pair((V) (((size_t) oldv) & ~set_bit), false);
  }
};

// Each thread maintains the log for the lock it is currently in using
// this variable.  It will be empty if not inside a lock.
static thread_local Log lg;

// executes the thunk f with log newlg
template <typename F>
auto inline with_log(Log newlg, F f) {
  Log holdlg = lg;
  lg = newlg;
  if constexpr (std::is_void_v<std::invoke_result_t<F>>) {
    f();
    lg = holdlg;
  } else {
    auto r = f();
    lg = holdlg;
    return r;
  }
}

template <typename F>
auto inline with_empty_log(F f) {
  return with_log(Log(), f);
}

// Skips a chunk of code if finished by another helper on the log.
// Uses a log entry to mark that it is completed.
template <typename F>
bool skip_if_done_no_log(F f) {
  if (lg.is_empty()) {f(); return true;}
  log_entry* l = lg.next_entry(); // get the next log entry
  if (l->load() == nullptr) { // check that not already completed
    f();
    // mark as completed
    l->store((void*) 1, std::memory_order::memory_order_release);
    return true;
  }
  return false;
}

template <typename F>
bool skip_if_done(F f) {
  if (lg.is_empty()) {f(); return true;}
  log_entry* l = lg.next_entry(); // get the next log entry
  void* v = l->load();
  if (v == nullptr) { // check that not already completed
    // run code
    f();
    // after running the code store the log location and log pointer
    size_t x = (((size_t) lg.count) << 48) | ((size_t) lg.vals);
    l->store((void*) x, std::memory_order::memory_order_release);
    return true;
  }
  // If we skip the code retrieve log location and pointer, so we can
  // "fast forward" to where would be in log if we ran the code.
  size_t x = (size_t) v;
  lg.count = (int) (x >> 48);
  lg.vals = (log_array*) (x & ((1ul << 48) - 1));
  return false;
}

// runs read only code without the log, and commits the result to the log
template <typename V, typename Thunk>
static V read_only(Thunk f) {
  // run f with an empty log and then commit its result
  V r = with_log(Log(), [&] {return f();});
  return lg.commit_value_safe(r).first;
}

} // namespace internal
} // namespace flck
