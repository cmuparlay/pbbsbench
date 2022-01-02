// The inteface for delaunay triangulation
#include "common/geometry.h"
#include "parlay/primitives.h"

using coord = double;
using point = point2d<coord>;

triangles<point> delaunay(parlay::sequence<point>& P);
