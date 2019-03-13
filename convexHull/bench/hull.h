#include "geometry.h"
#include "sequence.h"

using indexT = unsigned int;
using coord = double;
using point = point2d<coord>;

pbbs::sequence<indexT> hull(pbbs::sequence<point> const &S);

