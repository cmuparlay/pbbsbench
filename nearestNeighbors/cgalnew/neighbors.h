#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>

#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/get_time.h"
#include "../octTree/oct_tree.h"

int algorithm_version = 3; // just to interface with the timing code

using Kernel = CGAL::Simple_cartesian<double>;
using Point_3 = Kernel::Point_3;

using Traits = CGAL::Search_traits_3<Kernel>;
using Neighbor_search = CGAL::Orthogonal_k_neighbor_search<Traits>;
using Tree = Neighbor_search::Tree;
using Point_with_distance = Neighbor_search::Point_with_transformed_distance;

template<class vtx>
struct sorter{
  using o_tree = oct_tree<vtx>;
  using box = typename o_tree::box;

  parlay::sequence<vtx*> z_sort(parlay::sequence<vtx*> v, box b, double Delta){ 
    using indexed_point = typename o_tree::indexed_point; 
    size_t n = v.size();
    parlay::sequence<indexed_point> points;
    points = parlay::sequence<indexed_point>(n);
    parlay::parallel_for(0, n, [&] (size_t i){
      size_t p1 = o_tree::interleave_bits(v[i]->pt, b.first, Delta);
      indexed_point i1 = std::make_pair(p1, v[i]);
      points[i] = i1; 
    });
    auto less = [&] (indexed_point a, indexed_point b){
      return a.first < b.first;
    };
    auto x = parlay::sort(points, less);
    parlay::sequence<vtx*> v3; 
    v3 = parlay::sequence<vtx*>(n);
    parlay::parallel_for(0, n, [&] (size_t i){
      v3[i] = x[i].second; 
    });
    return v3; 
  }

};

template <int max_k, class vtx>
void ANN(parlay::sequence<vtx*> &v, int k) {
  using o_tree = oct_tree<vtx>;
  using box = typename o_tree::box;
  timer t("ANN",true);
  
  const unsigned int N = v.size();

  // std::vector<Point_3> points;
  // points.reserve (N);
  // for (unsigned int i = 0; i < N; ++ i) {
  //   auto pt = v[i]->pt;
  //   points.push_back (Point_3(pt[0],pt[1],pt[2]));
  // }

  auto minmax = [&] (box x, box y) {
    return box(x.first.minCoords(y.first),
         x.second.maxCoords(y.second));};
  
  // uses a delayed sequence to avoid making a copy
  auto pts = parlay::delayed_seq<box>(N, [&] (size_t i) {
      return box(v[i]->pt, v[i]->pt);});
  box identity = pts[0];
  box b = parlay::reduce(pts, parlay::make_monoid(minmax,identity));

  double Delta = 0;
  for (int i = 0; i < N; i++) 
    Delta = std::max(Delta, b.second[i] - b.first[i]); 

  sorter<vtx> S;
  parlay::sequence<vtx*> v1 = S.z_sort(v, b, Delta);

  auto points = parlay::tabulate(N, [&] (size_t i) {    
				      auto pt = v1[i]->pt;
				      return Point_3(pt[0],pt[1],pt[2]);});
  t.next("convert data");

  // Build tree in parallel
  Tree tree(points.begin(), points.end());
  tree.build<CGAL::Parallel_tag>();
  t.next("build tree");
  Point_3 a(0,0,0);

  // Query tree in parallel
  auto neighbors = parlay::tabulate(N, [&] (size_t s) {
                         // Neighbor search can be instantiated from
                         // several threads at the same time
                         Neighbor_search search (tree, points[s], k);
			 parlay::sequence<Point_3> ngh;
                         ngh.reserve(k);

                         // neighbor search returns a set of pair of
                         // point and distance <Point_3,FT>, here we
                         // keep the points only
                         for (const Point_with_distance& pwd : search)
                           ngh.push_back (pwd.first);
			 return ngh;
				       });
  t.next("search tree");
}
