#include "integer_sort.h"

template <class E, class F>
int num_bits(pbbs::sequence<E> In, F g) {
  auto get_key = [&] (size_t i) {return g(In[i]);};
  auto keys = pbbs::delayed_seq<size_t>(In.size(), get_key);
  return pbbs::log2_up(pbbs::reduce(keys, pbbs::maxm<size_t>()));
}

template <class T>
pbbs::sequence<T> int_sort(pbbs::sequence<T> const &In, size_t bits) {
  pbbs::sequence<T> Out = pbbs::sequence<T>::no_init(In.size());
  auto f = [&] (T x) {return x;};
  if (bits == 0) bits = num_bits(In,f);
  pbbs::seq_radix_sort(In.slice(), Out.slice(), In.slice(), f, bits);
  return Out;
}

template <class E, class F>
pbbs::sequence<std::pair<E,F>>
int_sort(pbbs::sequence<std::pair<E,F>> const &In, size_t bits) {
  using P = std::pair<E,F>;
  pbbs::sequence<P> Out = pbbs::sequence<P>::no_init(In.size());
  auto f = [&] (P x) {return x.first;};
  if (bits == 0) bits = num_bits(In,f);
  pbbs::seq_radix_sort(In.slice(), Out.slice(), In.slice(), f, bits);
  return Out;
}
