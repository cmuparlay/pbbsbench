/*
	Based on papers:
	PAM: Parallel Augmented maps
	Yihan Sun, Daniel Ferizovic and Guy Blelloch
	PPoPP 2018
	
	Parallel Range, Segment and Rectangle Queries with Augmented Maps
	Yihan Sun and Guy Blelloch
	arXiv:1803.08621
*/

#include <algorithm>
#include <iostream>
#include <vector>
#include <climits>
#include <cstdlib>
#include "parlay/primitives.h"
#include "parlay/internal/get_time.h"
#include "pam/pam.h"
#include "sweep.h"
#include "range.h"

using namespace std;

struct RangeQuery {
  struct set_t {
    using key_t = coord;
    static bool comp(key_t a, key_t b) { return a < b;}
  };
  using coord_set = pam_set<set_t>;

  parlay::sequence<coord_set> ts;
  parlay::sequence<coord> xs;

  RangeQuery(parlay::sequence<point> const &points) {
    auto less = [] (point a, point b) {return a.x < b.x;};
    auto A = parlay::sort(points, less);
    
    xs = parlay::map(A, [] (point p) {return p.x;});
    auto ys = parlay::map(A, [] (point p) {return p.y;});
    
    auto insert = [&] (coord_set m, coord a) {
      return coord_set::insert(m, a);};

    auto build = [&] (coord* s, coord* e) {
      return coord_set(s,e);};

    auto fold = [&] (coord_set m1, coord_set m2) {
      return coord_set::map_union(m1, std::move(m2));};

    ts = sweep<coord_set>(ys, coord_set(), insert, build, fold);
  }

  long count_in_range(query q) {
    long l = std::lower_bound(xs.begin(), xs.end(), q.x1) - xs.begin();
    long r = std::lower_bound(xs.begin(), xs.end(), q.x2) - xs.begin();
    long left = ts[l].rank(q.y2) - ts[l].rank(q.y1);
    long right = ts[r].rank(q.y2) - ts[r].rank(q.y1);
    return right-left;
  }

  void clear() {
    ts.clear();
    coord_set::GC::finish();
  }
};

long range(Points const &points, Queries const &queries) {
  parlay::internal::timer t("range");
  RangeQuery r(points);
  t.next("build");
  long total = parlay::reduce(parlay::map(queries, [&] (query q) {
  	          return (long) r.count_in_range(q);}));
  t.next("query");
  r.clear();
  t.next("clear");
  return total;
}
