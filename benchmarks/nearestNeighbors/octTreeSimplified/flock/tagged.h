#pragma once
#include <parlay/parallel.h>
#include <atomic>

// *****************************
// The following is used to tag values to avoid the ABA problem
// It uses the top 16 bits of a 64 bit value.
// It properly handles wrap around so tags are never reused if
// a thread has a pointer to them.
// It uses an announcement array to indicate that you have a pointer
// to a tag.
// tagged<V> adds a tag to type V.
// For a type V, and its tagged type TV, it supports
//    init(V) -> TV : add a tag
//    value(TV) -> V : strip off the tag
//    cas(std::atomic<TV*>, TV, V) -> bool :
//      compare and swap an untagged value V into the atomic if
//      it contains the tagged value TV.  return if successful
// The CAS is ensured to be ABA free.
// *****************************

namespace flck {
  namespace internal {
    
// The announcement array with one announcement per thread
struct write_annoucements {
  std::atomic<size_t>* announcements;
  const int stride = {16};
  write_annoucements() {
    announcements = new std::atomic<size_t>[parlay::num_workers()*stride];
  }
  ~write_annoucements() {
    delete[] announcements;
  }
  std::vector<size_t> scan() {
    std::vector<size_t> announced_tags;
    for(int i = 0; i < parlay::num_workers()*stride; i += stride)
      announced_tags.push_back(announcements[i]);
    return announced_tags;
  }
  void set(size_t val) {
    int id = parlay::worker_id();
    announcements[id*stride] = val;}
  void clear() {
    int id = parlay::worker_id();
    announcements[id*stride].store(0, std::memory_order::memory_order_release);}
};

write_annoucements announce_write = {};

// A wrapper to tag a value (either a pointer or a value with up to 48 bits).
template <typename V>
struct tagged {
private:
  using IT = size_t;
  static constexpr int tag_bits = 16; // number of bits to use for tag (including panic bit)
  static constexpr IT top_bit = (1ul << 63);
  static constexpr IT cnt_bit = (1ul << (64-tag_bits+1));
  static constexpr IT panic_bit = (1ul << (64-tag_bits));
  static constexpr IT data_mask = panic_bit - 1;
  static constexpr IT cnt_mask = ~data_mask;
  static inline IT add_tag(IT oldv, IT newv) {
    return newv | (oldv & cnt_mask);
  }
  static inline IT next(IT oldv, V newv) { // old version, doesn't handle overflow
    return ((IT) newv) | inc_tag(oldv); 
  }
  static inline IT next(IT oldv, V newv, IT addr) {
    // return next(oldv, newv);
    IT new_count = inc_tag(oldv);

    bool panic = false;
    if((oldv & top_bit) != (new_count & top_bit) || // overflow, unlikely
       (oldv & panic_bit) != 0) { // panic bit set, unlikely
      // if((oldv & top_bit) != (new_count & top_bit)) std::cout << "overflow\n";
      // if((oldv & panic_bit) != 0) std::cout << "panic bit set\n";
      for(IT ann : announce_write.scan()) { // check if we have to panic
        if((ann & data_mask) == (addr & data_mask) && // same mutable_val obj
           (ann & top_bit) == (new_count & top_bit) && // same half of the key range
           (ann & cnt_mask) >= (new_count & cnt_mask & ~panic_bit)) { 
          panic = true;
          break;
        }
      }
    }

    if(panic) { // unlikely
      std::vector<IT> announced_tags = announce_write.scan();
      while(true) { // loop until new_count is not announced
        bool announced = false;
        for(IT ann : announced_tags) {
          if((ann & data_mask) == (addr & data_mask) &&  // same mutable_val obj
             (ann & cnt_mask) == new_count) {
            announced = true;
            break;
          }
        }
        if(!announced) return ((IT) newv) | (new_count | panic_bit);
        new_count = inc_tag(new_count);
      }
    } else return ((IT) newv) | (new_count & ~panic_bit);
  }
  static inline IT inc_tag(IT oldv) {
    IT new_count = (oldv & cnt_mask) + cnt_bit;
    return ((new_count == 0) ? cnt_bit : new_count); // avoid using 0
  }

  // requires newV is already tagged
  static bool cas_tagged_(std::atomic<IT> &loc, IT oldv, IT newv, bool aba_free=false) {
    if (lg.is_empty() || aba_free) {
      return loc.compare_exchange_strong(oldv, newv);
    } else {
      bool r = false;
      // announce the location and tag been written
      announce_write.set(add_tag(oldv, (IT) &loc));
      skip_if_done_no_log([&] { // skip both for correctness, and efficiency
	  r = loc.compare_exchange_strong(oldv, newv);});
      // unannounce the location
      announce_write.clear();
      return r;
    }
  }
public:
  static inline IT init(V v) {return cnt_bit | (IT) v;}
  static inline V value(IT v) {return (V) (v & data_mask);}
  static inline IT get_tag(IT v) {return v & cnt_mask;}

  // a safe cas that assigns the new value a tag that no concurrent cas
  // on the same location has in its old value
  static bool cas(std::atomic<IT> &loc, IT oldv, V v, bool aba_free=false) {
    IT newv = next(oldv, v, (IT) &loc);
    return cas_tagged_(loc, oldv, newv, aba_free);
  }

  // a safe cas that assigns the new value a tag that no concurrent cas
  // on the same location has in its old value
  static bool cas_with_same_tag(std::atomic<IT> &loc, IT oldv, V v, bool aba_free=false) {
    IT newv = add_tag(oldv, (IT) v);
    return cas_tagged_(loc, oldv, newv, aba_free);
  }

};
  } // namespace internal
} // namespace flck
