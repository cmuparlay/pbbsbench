#include "parlay/primitives.h"

using charseq = parlay::sequence<char>;

charseq build_index(charseq const &s, charseq const &doc_start,
		    bool verbose);


