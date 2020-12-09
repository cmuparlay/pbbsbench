#include "parlay/primitives.h"

using charseq = parlay::sequence<unsigned char>;
using result_type = std::tuple<size_t,size_t,size_t>;
  
result_type lrs(charseq const &s);

