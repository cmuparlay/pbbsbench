#pragma once

#include <parlay/parallel.h> // needed for worker_id

#include <atomic>
#include<thread>
#include<optional>

thread_local bool in_lock;
struct lock {
private:
  flck::lock lck;
public:
  bool is_locked() { return lck.is_locked();}
  bool is_self_locked() { return lck.is_self_locked();}
  
  template <typename Thunk>
  auto try_lock(Thunk f) {
#ifdef Immediate
    if (current_transaction != nullptr) {
      bool taken = false;
      auto r = lck.try_lock_result(f, &taken); // do not release lock
      if (taken) current_transaction->lock_log.add(&lck);
      return r.has_value() && r.value();
    } else return lck.try_lock(f);
#else
    if (current_transaction != nullptr) {
      //if (lck.is_locked() && !lck.is_self_locked()) return false;
      current_transaction->lock_log.add(&lck);
      bool tmp = in_lock;
      in_lock = true;
      bool result = f();
      in_lock = tmp;
      return result;
    } else return lck.try_lock(f);
#endif
  }

  //   template <typename T, typename E>
  //   bool validate(T* v, E expected) {
  //     if (current_transaction != nullptr) {
  // #ifdef Validate
  //       if (try_count < 20) // anywhere from 10-100 seems good
  // 	return v->validate(expected);
  //       else // if too many tries give up on validation and take locks
  // 	return try_lock([=] {return (v->load() == expected);});
  // #else
  //       return try_lock([=] {return (v->load() == expected);});
  // #endif
  //     } else return true;
  //   }

  template <typename Thunk>
  auto try_lock_no_delay(Thunk f) {return try_lock(f);}
};

template <typename F>
bool run_with_locks_rec(int i, trans_descriptor* descriptor, F f) {
  if (i == descriptor->lock_log.size()) return f();
  return descriptor->lock_log[i]->try_lock([&] {
					     return run_with_locks_rec(i+1, descriptor, f);});
}

template <typename F>
bool run_with_logged_locks(trans_descriptor* descriptor, F f) {
  // sorting actually slows things down slightly (5%) under high contention
  std::sort(descriptor->lock_log.data(), descriptor->lock_log.data() + descriptor->lock_log.size());
  // long max_count = lock_log.vals[0]->count();
  // int j=0;
  // for (int i=1; i < lock_log.cnt; i++) {
  //   long cnt = lock_log.vals[i]->count();
  //   if (cnt > max_count) {j = i; max_count = cnt;}
  // }
  // std::swap(lock_log.vals[j], lock_log.vals[lock_log.cnt-1]);
  // if (false)
  //   for (int i=0; i < validate_ptr_log.cnt; i++) {
  // 	auto [location, expected] = validate_ptr_log.vals[i];
  // 	if (!validate_pointer(location, expected)) return false;
  //   }
  if (false)
    for (int i=0; i < descriptor->lock_log.size(); i++) {
      if (descriptor->lock_log[i]->is_locked() &&
	  !descriptor->lock_log[i]->is_self_locked())
	return false;
    }
  return run_with_locks_rec(0, descriptor, f);
}
