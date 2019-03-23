#define NOTMAIN 1
#include "strings/suffix_array.h"
#include "SA.h"

pbbs::sequence<indexT> suffixArray(pbbs::sequence<unsigned char> const &s) {
  return pbbs::suffix_array<indexT>(s);
}
