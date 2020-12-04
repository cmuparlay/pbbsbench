#include "parlay/primitives.h"

using charseq = parlay::sequence<char>;
using result_type = std::pair<charseq,size_t>;
parlay::sequence<result_type> wordCounts(charseq const &s, bool verbose);


