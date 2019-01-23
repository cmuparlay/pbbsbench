#include "pbbslib/integer_sort.h"

template <class E>
void int_sort(E* A, unsigned int n) {
  sequence<E> In(A,n);
  auto f = [&] (E x) {return x;};
  pbbs::integer_sort(In, f);
}

template <class E, class F>
void int_sort(std::pair<E,F>* A, unsigned int n) {
  sequence<std::pair<E,F>> In(A,n);
  auto f = [&] (std::pair<E,F> x) {return x.first;};
  pbbs::integer_sort(In, f);
}
