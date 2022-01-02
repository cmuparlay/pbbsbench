#include "parlay/sequence.h"
#include "parlay/internal/get_time.h"
#include "algorithm/suffix_array.h"
#include "algorithm/lcp.h"

using charseq = parlay::sequence<unsigned char>;
using result_type = std::tuple<size_t,size_t,size_t>;

// returns
//  1) the length of the longest match
//  2) start of the first string in s
//  3) start of the second string in s
template <typename IntType>
result_type lrs_(charseq const &s) {
  parlay::internal::timer t("lrs", true);

  parlay::sequence<IntType> sa = suffix_array<IntType>(s);
  t.next("suffix array");

  parlay::sequence<IntType> lcps = lcp(s, sa);
  t.next("lcps");

  size_t idx = parlay::max_element(lcps, std::less<IntType>())-lcps.begin();
  t.next("max element");
    
  return result_type(lcps[idx],sa[idx],sa[idx+1]);
}

result_type lrs(charseq const &s) {
  return lrs_<unsigned int>(s);
}
