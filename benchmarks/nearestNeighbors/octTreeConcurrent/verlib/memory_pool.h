namespace verlib {
  template <typename xT>
  struct memory_pool {
    using T = xT;
    flck::memory_pool<T> pool;
    void reserve(size_t n) { pool.reserve(n);}
    void shuffle(size_t n) { pool.shuffle(n);}
    void stats() { pool.stats();}
    void clear() { pool.clear();}
    void destruct(T* ptr) { pool.destruct(ptr); }

    bool is_helping() {
#ifdef NoHelp
      return false;
#else
      return flck::internal::helping;
#endif
    }
  
    void retire_on_abort(trans_descriptor* descriptor, T* ptr) {
#ifdef TransUndo
      if (descriptor != nullptr && !is_helping()) {
	bool* x = pool.retire(ptr);
	flck::internal::undo_retire(x);
	descriptor->allocated_log.add(x);
      }
#endif
    }

    void retire_on_success(trans_descriptor* descriptor, T* ptr) {
#ifdef TransUndo
      bool* x = pool.retire(ptr);
      if (descriptor != nullptr && !is_helping())
	descriptor->retired_log.add(x);
#else
      pool.retire(ptr);
#endif
    }

    template <typename ... Args>
    T* new_obj(Args... args) {
      T* ptr = pool.new_obj(args...);
      retire_on_abort(current_transaction, ptr);
      return ptr;
    }

    template <typename F, typename ... Args>
    T* new_init(F f, Args... args) {
      T* x = new_obj(args...);
      f(x);
      return x;
    }

    void retire(T* ptr) {
      retire_on_success(current_transaction, ptr);
    }
  };
}
