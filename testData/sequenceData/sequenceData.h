// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2010 Guy Blelloch and the PBBS team
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

#include <iostream>
#include <fstream>
#include <math.h>
#include "parlay/primitives.h"
#include "parlay/random.h"

namespace dataGen {

  struct payload {
    double key;
    double payload[2];
  };

  class payloadCmp : public std::binary_function <payload, payload, bool> {
  public:
    bool operator()(payload const& A, payload const& B) const {
      return A.key<B.key;
    }
  };

  template<typename T>
  parlay::sequence<T> randIntRange(size_t s, size_t e, size_t m) {
    parlay::random r(0);
    auto In = parlay::tabulate(e-s, [&] (size_t i) -> T {return
	  r.ith_rand(i+s)%m;});
    return In;
  }

  template<typename T>
  parlay::sequence<T> rand(size_t s, size_t e) {
    parlay::random r(0);
    auto In = parlay::tabulate(e-s, [&] (size_t i) -> T{return
	  r.ith_rand(i+s);});
    return In;
  }

  template <class T>
  parlay::sequence<T> almostSorted(size_t s, size_t e, size_t swaps) { 
    size_t n = e - s;
    parlay::random r(0);
    auto A = parlay::tabulate(n, [&] (long i) -> T {return i;});
    for (size_t i = s; i < s+swaps; i++)
      std::swap(A[r.ith_rand(2*i)%n],A[r.ith_rand(2*i+1)%n]);
    return A;
  }

  template <class T>
  parlay::sequence<T> same(size_t n, T v) {
    parlay::sequence<T> A(n,v);
    return A;
  }

  template <class T>
  parlay::sequence<T> expDist(size_t s, size_t e) { 
    size_t n = e - s;
    size_t lg = parlay::log2_up(n)+1;
    parlay::random r(0);
    auto A = parlay::tabulate(n, [&] (long i) -> T {
      size_t range = (1 << (r.ith_rand(2*(i+s))%lg));
      return r.ith_rand((size_t)(range + r.ith_rand(2*(i+s)+1)%range))%n;
      });
    return A;
  }

};
