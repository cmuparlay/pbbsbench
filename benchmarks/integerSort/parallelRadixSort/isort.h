#include "parlay/primitives.h"

template <class T>
auto int_sort(parlay::slice<T*,T*> In, size_t bits) {
  auto f = [&] (T x) {return x;};
  return parlay::internal::integer_sort(parlay::make_slice(In), f, bits);
}

template <class E, class F>
auto int_sort(parlay::slice<std::pair<E,F>*, std::pair<E,F>*> In, size_t bits) {
  auto f = [&] (std::pair<E,F> x) {return x.first;};
  return parlay::internal::integer_sort(parlay::make_slice(In), f, bits);
}
