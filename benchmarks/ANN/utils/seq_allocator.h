// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef FINE_ALLOCATOR
#define FINE_ALLOCATOR

#include "parlay/alloc.h"
#include <cmath>

// these are bucket sizes used by the fine allocator.
inline std::vector<size_t> default_sizes() {
  // size_t log_min_size = 4;
  // size_t log_max_size = parlay::log2_up(getMemorySize()/64);

  std::vector<size_t> sizes;
  size_t rounded_current_size = 8;
  double current_size = 8;
  while(rounded_current_size < 300) {
    if (((size_t)current_size) > rounded_current_size) {
      rounded_current_size = current_size;
      sizes.push_back(rounded_current_size);
    }
  	current_size = current_size*1.1;
  }
  return sizes;
}


namespace fine{

	extern inline parlay::pool_allocator& get_fine_default_allocator() {
	  static parlay::pool_allocator fine_default_allocator(default_sizes());
	  return fine_default_allocator;
	}

	// pair of total currently used space, and total unused space the allocator has in reserve
	extern inline std::pair<size_t,size_t> memory_usage() {
	  return get_fine_default_allocator().stats();
	}

	// pair of total currently used space, and total unused space the allocator has in reserve
	extern inline void memory_clear() {
	  return get_fine_default_allocator().clear();
	}
} //namespace fine


// ****************************************
// Following Matches the c++ Allocator specification (minimally)
// https://en.cppreference.com/w/cpp/named_req/Allocator
// Can therefore be used for containers, e.g.:
//    std::vector<int, parlay::allocator<int>>
// ****************************************

template <typename T>
struct fine_allocator {
  using value_type = T;
  T* allocate(size_t n) {
    return (T*) fine::get_fine_default_allocator().allocate(n * sizeof(T));
  }
  void deallocate(T* ptr, size_t n) {
    fine::get_fine_default_allocator().deallocate((void*) ptr, n * sizeof(T));
  }

  constexpr fine_allocator() = default;
  template <class U> constexpr fine_allocator(const fine_allocator<U>&) noexcept { }
};

template <class T, class U>
bool operator==(const fine_allocator<T>&, const fine_allocator<U>&) { return true; }
template <class T, class U>
bool operator!=(const fine_allocator<T>&, const fine_allocator<U>&) { return false; }

constexpr size_t size_offset = 1; // in size_t sized words

// needs to be at least size_offset * sizeof(size_t)
inline size_t header_size(size_t n) { // in bytes
  return (n >= 1024) ? 64 : (n & 15) ? 8 : (n & 63) ? 16 : 64;
}

// allocates and tags with a header (8, 16 or 64 bytes) that contains the size
extern inline void* p_malloc(size_t n) {
  size_t hsize = header_size(n);
  void* ptr = fine::get_fine_default_allocator().allocate(n + hsize);
  void* r = (void*) (((char*) ptr) + hsize);
  *(((size_t*) r) - size_offset) = n; // puts size in header
  return r;
}

// reads the size, offsets the header and frees
extern inline void p_free(void* ptr) {
  size_t n = *(((size_t*) ptr) - size_offset);
  size_t hsize = header_size(n);
  if (hsize > (1ull << 48)) {
    std::cout << "corrupted header in my_free" << std::endl;
    throw std::bad_alloc();
  }
  fine::get_fine_default_allocator().deallocate((void*) (((char*) ptr) - hsize),
                                               n + hsize);
}

#endif