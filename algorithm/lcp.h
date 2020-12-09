#include "../parlay/sequence.h"
#include "../parlay/internal/get_time.h"
#include "range_min.h"

//  The suffix array SA are indices into the string s
template <class Seq1, class Seq2>
auto lcp(Seq1 const &s_, Seq2 const &SA_)
  -> parlay::sequence<typename Seq2::value_type>
{
  parlay::internal::timer t("LCP", false);
  auto s = parlay::make_slice(s_); 
  auto SA = parlay::make_slice(SA_);
  using Uint = typename Seq2::value_type;
  size_t len = 111;
  size_t n = SA.size();
  t.next("init");
    
  // compare first len characters of adjacent strings from SA.
  auto L_ = parlay::tabulate(n-1, [&] (size_t i) -> Uint {
	       size_t j = 0;
	       size_t max_j = std::min(len, n - SA[i]);
	       while (j < max_j && (s[SA[i]+j] == s[SA[i+1]+j])) j++;
	       return (j < len) ? j : n;
	    });
  auto L = parlay::make_slice(L_);
  t.next("head");

  // keep indices for which we do not yet know their LCP (i.e. LCP >= len)
  auto remain = parlay::pack_index(parlay::map(L, [&] (Uint l) {return l == n;}));
  t.next("pack");

  if (remain.size() == 0) return L_;

  // an inverse permutation for SA
  parlay::sequence<Uint> ISA_(n);
  auto ISA = parlay::make_slice(ISA_);
  parlay::parallel_for(0, n, [&] (size_t i) {
			       ISA[SA[i]] = i;});
  t.next("inverse");

  // repeatedly double len determining LCP by joining next len chars
  // invariant: before i^th round L contains correct LCPs less than len
  //            and n for the rest of them
  //            remain holds indices of the rest of them (i.e., LCP[i] >= len)
  //      after round, len = 2*len and invariant holds for the new len
  do {
    auto rq = make_range_min(L, std::less<Uint>(), 111);
    t.next("make range");

    // see if next len chars resolves LCP
    // set L for those that do, and keep those that do not for next round
    remain = parlay::filter(remain, [&] (Uint i) {
		if (SA[i] + len >= n) {L[i] = len; return false;};
		Uint i1 = ISA[SA[i]+len];
		Uint i2 = ISA[SA[i+1]+len];
		size_t l = L[rq.query(i1, i2-1)];
		if (l < len) {L[i] = len + l; return false;}
		else return true;
	     });
    t.next("filter");
    len *= 2;
  } while (remain.size() > 0);

  return L_;
}


