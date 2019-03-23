#define NOTMAIN 1
#include "strings/suffix_array.h"
#include "SA.h"

pbbs::sequence<unsigned int> suffixArray(pbbs::sequence<unsigned char> const &s) {
  return pbbs::suffix_array<unsigned int>(s);
}
