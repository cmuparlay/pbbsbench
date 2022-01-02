#include "common/geometry.h"
#include "parlay/sequence.h"

using coord = double;
using point = point3d<double>;
using index_t = int; // if more than 2 billion rays, then should be long

parlay::sequence<index_t> rayCast(triangles<point> const &Tri,
				  parlay::sequence<ray<point>> const &rays, bool verbose);
