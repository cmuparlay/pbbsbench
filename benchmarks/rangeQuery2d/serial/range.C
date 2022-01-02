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
#include <limits>
#include <cfloat>
#include "parlay/primitives.h"
#include "parlay/internal/get_time.h"
#include "pam/pam.h"
#include "sweep.h"
#include "range.h"

struct RangeQuery {
  struct set_t {
    using key_t = point;
    static bool comp(key_t a, key_t b) { 
	   if ((a.y < b.y) || ( (a.y == b.y) && (a.x < b.x) ) ) return true; else return false;
	}
  };
  using coord_set = pam_set<set_t>;

  parlay::sequence<coord_set> ts;
  parlay::sequence<coord> xs;

  RangeQuery(parlay::sequence<point> const &points) {
    auto less = [] (point a, point b) {return a.x < b.x;};
    auto A = parlay::sort(points, less);
    
    xs = parlay::map(A, [] (point p) {return p.x;});
    //auto ys = parlay::map(A, [] (point p) {return p;});
	parlay::sequence<point> ys = points;
    
    auto insert = [&] (coord_set m, point a) {
      return coord_set::insert(m, a);};

    auto build = [&] (point* s, point* e) {
      return coord_set(s,e);};

    auto fold = [&] (coord_set m1, coord_set m2) {
      return coord_set::map_union(m1, std::move(m2));};

    ts = sweep<coord_set>(ys, coord_set(), insert, build, fold, 1);
	
  }

  long count_in_range(query q) {
    long l = std::lower_bound(xs.begin(), xs.end(), q.x1) - xs.begin();
    long r = std::lower_bound(xs.begin(), xs.end(), q.x2) - xs.begin();
    long left = ts[l].rank(point(DBL_MAX, q.y2)) - ts[l].rank(point(-DBL_MAX, q.y1));
    long right = ts[r].rank(point(DBL_MAX, q.y2)) - ts[r].rank(point(-DBL_MAX, q.y1));
    return right-left;
  }

  void clear() {
    ts.clear();
    coord_set::GC::finish();
  }
};

long range(Points const &points, Queries const &queries, bool verbose) {
  parlay::internal::timer t("range", verbose);
  RangeQuery r(points);
  t.next("build");
  long total = parlay::reduce(parlay::map(queries, [&] (query q) {
  	          return (long) r.count_in_range(q);}));
  /*auto result_s = parlay::map(queries, [&] (query q) {
  	          return (long) r.count_in_range(q);});
  long total = parlay::reduce(result_s);*/
  t.next("query");

#ifdef CHECK
  //check
  //int num_queries = queries.size();
  int num_queries = 10; //only check 10 of them 
  int n = points.size();
  size_t total_check = 0;
	for (int i = 0; i < num_queries; i++) {
		//cout << "query: " << queries[i].x1 << " " << queries[i].y1 << " " << queries[i].x2 << " " << queries[i].y2 << endl;
		size_t t = 0;
		for (int j = 0; j < n; j++) {
			 if ((points[j].x > queries[i].x1) && (points[j].x < queries[i].x2) && (points[j].y > queries[i].y1) && (points[j].y < queries[i].y2)) t++;
		}
		long res = r.count_in_range(queries[i]);
		cout << "query " << i << ", naive: " << t << ", PAM: " << res << endl;
		total_check += t;
	}
  cout << "total: " << total << endl;
  cout << "total_check: " << total_check << endl;
#endif

  r.clear();
  t.next("clear");
  return total;
}
