#include "pbbslib/integer_sort.h"

template <class T>
pbbs::sequence<T> int_sort(pbbs::sequence<T> const &In, size_t bits) {
  auto f = [&] (T x) {return x;};
  return pbbs::integer_sort(In, f, bits);
}

template <class E, class F>
pbbs::sequence<std::pair<E,F>> int_sort(pbbs::sequence<std::pair<E,F>> const &In, size_t bits) {
  auto f = [&] (std::pair<E,F> x) {return x.first;};
  return pbbs::integer_sort(In, f, bits);
}
