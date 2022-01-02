#include "common/geometry.h"
#include "parlay/primitives.h"

using coord = double;
using point = point3d<coord>;
using vect = vector3d<coord>;

class particle {
public:
  point pt;
  vect force;
  double mass;
  particle(point p, double m) : pt(p), mass(m) {}
};

void nbody(parlay::sequence<particle*> &particles);
