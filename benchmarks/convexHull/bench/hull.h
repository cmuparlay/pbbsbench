#include "common/geometry.h"
#include "parlay/primitives.h"

using indexT = unsigned int;
using coord = double;
using point = point2d<coord>;

parlay::sequence<indexT> hull(parlay::sequence<point> const &S);

