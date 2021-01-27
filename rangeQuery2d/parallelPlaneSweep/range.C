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

  using entry_t = pair<coord, weight>;
  
  struct map_t {
    using key_t = coord;
    using val_t = weight;
    static bool comp(key_t a, key_t b) { return a < b;}
    using aug_t = weight;
    static aug_t get_empty() {return 0;}
    static aug_t from_entry(key_t k, val_t v) {return v;}
    static aug_t combine(aug_t a, aug_t b) {return a + b;}
  };

  using c_map = aug_map<map_t>;
  parlay::sequence<c_map> ts;
  parlay::sequence<coord> xs;
  size_t n;

  RangeQuery(parlay::sequence<point> const &points) {
    n = points.size();
    c_map::reserve(36 * n);

    auto less = [] (point a, point b) {return a.x < b.x;};
    auto A = parlay::sort(points, less);
    
    xs = parlay::sequence<coord>::uninitialized(n);
    auto ys = parlay::sequence<entry_t>::uninitialized(n);
    parlay::parallel_for (0, n, [&] (size_t i) {
      xs[i] = A[i].x;
      ys[i] = entry_t(A[i].y, A[i].w);
      });

    auto plus = [] (weight a, weight b) {return a + b;};
    
    auto insert = [&] (c_map m, entry_t a) {
      return c_map::insert(m, a, plus);
    };

    auto build = [&] (entry_t* s, entry_t* e) {
      return c_map(s,e,plus);
    };

    auto fold = [&] (c_map m1, c_map m2) {
      return c_map::map_union(m1, std::move(m2), plus);};

    ts = sweep<c_map>(ys, c_map(), insert, build, fold);
  }

  int get_index(const coord q) {
    int l = 0, r = n;
    int mid = (l+r)/2;
    if (xs[0]>q) return -1;
    while (l<r-1) {
      if (xs[mid] == q) break;
      if (xs[mid] < q) l = mid;
      else r = mid;
      mid = (l+r)/2;
    }
    return mid;
  }

  weight count_in_range(query q) {
    int l = get_index(q.x1);
    int r = get_index(q.x2);
    weight left = (l<0) ? 0.0 : ts[l].aug_range(q.y1, q.y2);
    weight right = (r<0) ? 0.0 : ts[r].aug_range(q.y1, q.y2);
    return right-left;
  }

  static void print_allocation_stats() {
    cout << "allocation stats:" ;  c_map::GC::print_stats();
  }

  static void finish() {
    c_map::GC::finish();
  }
};

long range(Points const &points, Queries const &queries) {
  parlay::internal::timer t("range");
  RangeQuery r(points);
  t.next("build");
  query q = queries[0];
  long total = parlay::reduce(parlay::map(queries, [&] (query q) {
  	          return (long) r.count_in_range(q);}));
  t.next("query");
  return total;
}
