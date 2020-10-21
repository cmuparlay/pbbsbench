#include "parlay/primitives.h"
#include "parlay/parallel.h"

template <class Seq, class F>
int num_bits(Seq In, F g) {
  auto get_key = [&] (size_t i) {return g(In[i]);};
  auto keys = parlay::delayed_seq<size_t>(In.size(), get_key);
  return parlay::log2_up(parlay::reduce(keys, parlay::maxm<size_t>()));
}

template <class T>
auto int_sort(parlay::slice<T*,T*> In, size_t bits) {
  auto Out = parlay::sequence<T>::uninitialized(In.size());
  auto f = [&] (T x) {return x;};
  if (bits == 0) bits = num_bits(In,f);
  parlay::internal::seq_radix_sort<std::true_type,parlay::uninitialized_relocate_tag>
       (In, make_slice(Out), In, f, bits);
  return Out;
}

template <class E, class F>
auto int_sort(parlay::slice<std::pair<E,F>*, std::pair<E,F>*> In, size_t bits) {
  using P = std::pair<E,F>;
  auto Out = parlay::sequence<P>::uninitialized(In.size());
  auto f = [&] (P x) {return x.first;};
  if (bits == 0) bits = num_bits(In,f);
  parlay::internal::seq_radix_sort<std::true_type,parlay::uninitialized_relocate_tag>
       (In, make_slice(Out), In, f, bits);
  return Out;
}
