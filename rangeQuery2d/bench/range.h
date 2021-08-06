using coord = double;
struct point {
  coord x, y;
  point(coord x, coord y) : x(x), y(y) {}
  point() {}
};
struct query {coord x1, x2, y1, y2;};
using Points = parlay::sequence<point>;
using Queries = parlay::sequence<query>;

long range(Points const &points, Queries const &queries, bool verbose);
