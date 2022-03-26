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
  size_t current_size = 4;
  while(current_size < 300){
  	current_size = (size_t) current_size*1.2;
  	sizes.push_back(current_size);
  }
  return sizes;
}



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
    return (T*) get_fine_default_allocator().allocate(n * sizeof(T));
  }
  void deallocate(T* ptr, size_t n) {
    get_fine_default_allocator().deallocate((void*) ptr, n * sizeof(T));
  }

  constexpr fine_allocator() = default;
  template <class U> constexpr fine_allocator(const allocator<U>&) noexcept { }
};

#endif