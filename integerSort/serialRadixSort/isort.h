#include "pbbslib/integer_sort.h"

template <class E, class F>
int num_bits(sequence<E> In, F g) {
  auto max = [&] (size_t a, size_t b) {return std::max(a,b);};
  auto get_key = [&] (size_t i) {return g(In[i]);};
  auto keys = make_sequence<size_t>(In.size(), get_key);
  return pbbs::log2_up(pbbs::reduce(keys, max));
}

template <class E>
void int_sort(E* A, unsigned int n) {
  sequence<E> In(A,n);
  sequence<E> Out = sequence<E>::alloc_no_init(n);
  auto f = [&] (E x) {return x;};
  int bits = num_bits(In,f);
  pbbs::seq_radix_sort(In, Out, f, bits, true);
}

template <class E, class F>
void int_sort(std::pair<E,F>* A, unsigned int n) {
  using P = std::pair<E,F>;
  sequence<P> In(A,n);
  sequence<P> Out = sequence<P>::alloc_no_init(n);
  auto f = [&] (P x) {return x.first;};
  int bits = num_bits(In,f);
  pbbs::seq_radix_sort(In, Out, f, bits, true);
}
