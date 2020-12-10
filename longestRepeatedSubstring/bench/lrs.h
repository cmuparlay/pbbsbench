#include "parlay/primitives.h"

using charseq = parlay::sequence<unsigned char>;

// returns
//  1) the length of the longest match
//  2) start of the first string in s
//  3) start of the second string in s
using result_type = std::tuple<size_t,size_t,size_t>;
  
result_type lrs(charseq const &s);

