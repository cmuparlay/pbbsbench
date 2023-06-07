// Currently not used

struct hold {
    atomic_val old_val;
    atomic_val new_val;
    trans_descriptor* descriptor;
    std::atomic<bool> acquired;
    constexpr static int bit_num = 63;
    static bool is_held(atomic_val x) {return (x >> bit_num) & 1;}
    static atomic_val tag_hold(hold* ptr) {
      return (TS) ptr | (1ul << bit_num);}
    static hold* untag_hold(atomic_val x) {
      return (hold*) (x & ~(1ul << bit_num)); }
    hold(atomic_val o, atomic_val n, trans_descriptor* d)
      : old_val(o), new_val(n), descriptor(d) {}
  };

  //static flck::memory_pool<hold> hold_pool;
  static flck::internal::acquired_pool<hold> hold_pool;

  void install_atomic_store(std::atomic<atomic_val>* loc,
			    atomic_val new_v,
			    trans_descriptor* descriptor) {
    atomic_val old_v = loc->load();
    hold* ptr = hold_pool.new_obj(old_v, (atomic_val) new_v, descriptor);
    *loc = hold::tag_hold(ptr);
  }

  void finalize_atomic_store(std::atomic<atomic_val>* loc, TS time_stamp) {
    assert(hold::is_held(*loc));
    auto hold_ptr = hold::untag_hold(*loc);
    assert(hold_ptr->descriptor->owner == flck::internal::current_id);
    assert(hold_ptr->descriptor->time_stamp == time_stamp);
    *loc = (time_stamp == failed) ? hold_ptr->old_val : hold_ptr->new_val;
    hold_pool.retire(hold_ptr);
  }

  V get_value(std::atomic<atomic_val>& v, bool validating=false) {
    atomic_val x = v.load();
    atomic_val y;
    hold* ptr;
    do {
      if (!hold::is_held(x)) return x;
      ptr = hold::untag_hold(x);
      hold_pool.acquire(ptr);
      y = x;
      x = v.load();
    } while (x != y);
    return resolve_descriptor(ptr->descriptor, ptr->old_val,
			      ptr->new_val, validating);
  }

  template <typename V> 
  struct atomic {
    std::atomic<atomic_val> v;
    atomic(V v) : v(v) {}
    std::optional<int> find_in_atomic_write_log() {
      for (int i=0; i < write_atomic_log.cnt; i++)
	if (write_atomic_log.vals[i].first == &v)
	  return std::optional<int>(i);
      return std::optional<int>();
    }
    V load() {
      if (current_transaction != nullptr) {
	std::optional<int> lid = find_in_atomic_write_log();
	if (lid.has_value())
	  return (V*) write_atomic_log.vals[*lid].second;
	if (in_lock) {
	  atomic_val cur = get_value(v);
	  validate_atomic_log.add(std::tuple(&v, cur));
	  return (V) cur;
	}
      }
      return (V) get_value(v);
    }
    void store(V new_v) {
      if (current_transaction == nullptr) v = new_v;
      else {
	std::optional<int> lid = find_in_atomic_write_log();
	if (lid.has_value())
	  write_atomic_log.vals[*lid].second = (atomic_val) new_v;
	else write_atomic_log.add(std::pair(&v, new_v));
      }
    }
    bool validate(V expected) { 
      std::optional<int> lid = find_in_atomic_write_log();
      if (lid.has_value()) return (write_atomic_log.vals[*lid].second == expected);
      atomic_val cur = get_value(v);
      if (cur != expected) return false; 
      validate_atomic_log.add(std::tuple(&v, expected));
      return true;
    }
    V operator=(V b) {store(b); return b; }
  };
