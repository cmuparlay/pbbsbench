using coord = double;
using weight = int;

struct point {
  coord x, y;
  weight w;
  point() {}
  point(coord x, coord y, weight w) : x(x), y(y), w(w) {}
};

struct query {
  coord x1, x2, y1, y2;
};

using Points = parlay::sequence<point>;
using Queries = parlay::sequence<query>;

long range(Points const &points, Queries const &queries);
