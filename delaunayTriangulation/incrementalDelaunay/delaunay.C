// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include <vector>
#include "parlay/primitives.h"
#include "parlay/random.h"
#include "common/geometry.h"
#include "common/get_time.h"
#include "common/topology.h"
#include "common/atomics.h"
#include "neighbors.h"
#include "delaunay.h"

using parlay::parallel_for;
using parlay::sequence;
using parlay::delayed_seq;
using parlay::tabulate;
using parlay::reduce;
using parlay::pack;
using parlay::make_monoid;
using parlay::random_permutation;
using parlay::internal::pack_out;
using std::cout;
using std::endl;

// if on verifies the Delaunay is correct 
#define CHECK 0

using vertex_t = vertex<point>;
using simplex_t = simplex<point>;
using triang_t = triangle<point>;
using vect = typename point::vector;

template <typename point>
struct Qs {
  vector<vertex<point>*> vertexQ;
  vector<simplex<point>> simplexQ;
  Qs() {
    vertexQ.reserve(50);
    simplexQ.reserve(50);
  }
};
using Qs_t = Qs<point>;

// *************************************************************
//    ROUTINES FOR FINDING AND INSERTING A NEW POINT
// *************************************************************

// Finds a vertex (p) in a mesh starting at any triangle (start)
// Requires that the mesh is properly connected and convex
simplex_t find(vertex_t *p, simplex_t start) {
  simplex_t t = start;
  while (1) {
    int i;
    for (i=0; i < 3; i++) {
      t = t.rotClockwise();
      if (t.outside(p)) {t = t.across(); break;}
    }
    if (i==3) return t;
    if (!t.valid()) return t;
  }
}

// Recursive routine for finding a cavity across an edge with
// respect to a vertex p.
// The simplex has orientation facing the direction it is entered.
//
//         a
//         | \ --> recursive call
//   p --> |T c 
// enter   | / --> recursive call
//         b
//
//  If p is in circumcircle of T then 
//     add T to simplexQ, c to vertexQ, and recurse
void findCavity(simplex_t t, vertex_t *p, Qs_t *q) {
  if (t.inCirc(p)) {
    q->simplexQ.push_back(t);
    t = t.rotClockwise();
    findCavity(t.across(), p, q);
    q->vertexQ.push_back(t.firstVertex());
    t = t.rotClockwise();
    findCavity(t.across(), p, q);
  }
}

// Finds the cavity for v and tries to reserve vertices on the 
// boundary (v must be inside of the simplex t)
// The boundary vertices are pushed onto q->vertexQ and
// simplices to be deleted on q->simplexQ (both initially empty)
// It makes no side effects to the mesh other than to X->reserve
void reserve_for_insert(vertex_t *v, simplex_t t, Qs_t *q) {
  // each iteration searches out from one edge of the triangle
  for (int i=0; i < 3; i++) {
    q->vertexQ.push_back(t.firstVertex());
    findCavity(t.across(), v, q);
    t = t.rotClockwise();
  }
  // the maximum id new vertex that tries to reserve a boundary vertex 
  // will have its id written.  reserve starts out as -1
  for (size_t i = 0; i < q->vertexQ.size(); i++) {
    //cout << "trying to reserve: " << (q->vertexQ)[i]->reserve << ", " << v->id << endl;
    pbbs::write_max(&((q->vertexQ)[i]->reserve), v->id, std::less<int>());
  }
}

// checks if v "won" on all adjacent vertices and inserts point if so
bool insert(vertex_t *v, simplex_t t, Qs_t *q) {
  bool flag = 0;
  for (size_t i = 0; i < q->vertexQ.size(); i++) {
    vertex_t* u = (q->vertexQ)[i];
    //cout << u->id << ", " << u->reserve << " : " << v->id << endl;
    if (u->reserve == v->id) u->reserve = -1; // reset to -1
    else flag = 1; // someone else with higher priority reserved u
  }
  if (!flag) {
    triang_t* t1 = v->t;  // the memory for the two new triangles
    triang_t* t2 = t1 + 1;  
    // the following 3 lines do all the side effects to the mesh.
    t.split(v, t1, t2);
    //cout << "just split: " << q->simplexQ.size() << endl;
    for (size_t i = 0; i<q->simplexQ.size(); i++) {
      //cout << "flipping: " << i << endl;
      (q->simplexQ)[i].flip();
    }
  }
  q->simplexQ.clear();
  q->vertexQ.clear();
  return flag;
}

// *************************************************************
//    CHECKING THE TRIANGULATION
// *************************************************************

void check_delaunay(sequence<triang_t> &Triangles, size_t boundary_size) {
  size_t n = Triangles.size();
  sequence<size_t> boundary_count(n, 0);
  parallel_for (0, n, [&] (size_t i) {
    if (Triangles[i].initialized >= 0) {
      simplex_t t = simplex(&Triangles[i], 0);
      for (int i=0; i < 3; i++) {
	simplex_t a = t.across();
	if (a.valid()) {
	  vertex_t* v = a.rotClockwise().firstVertex();
	  if (!t.outside(v)) {
	    cout << "Inside Out: "; v->pt.print(); t.print();}
	  if (t.inCirc(v)) {
	    cout << "In Circle Violation: "; v->pt.print(); t.print(); }
	} else boundary_count[i]++;
	t = t.rotClockwise();
      }
    } });
  if (boundary_size != reduce(boundary_count))
    cout << "Wrong boundary size: should be " << boundary_size 
	 << " is " << reduce(boundary_count) << endl;
}

// *************************************************************
//    CREATING A BOUNDING CIRCULAR REGION AND FILL WITH INITIAL SIMPLICES
// *************************************************************

// P is the set of points to bound and n the number
// boundary_size is the number of points to put on the boundary
// V is a sequence of vertices, which the new vertices are added to, at end
// T is a sequence of triangles, which the new triangles are added to, at end
// one of the triangles is returned as an ordered simplex
void generate_boundary(sequence<point> const &P,
		       size_t boundary_size,
		       sequence<vertex_t> &V,
		       sequence<triang_t> &T) {

  size_t n = P.size();
  auto min = [] (point x, point y) { return x.minCoords(y);};
  auto max = [] (point x, point y) { return x.maxCoords(y);};
  point identity = P[0];
  point min_corner = reduce(P, make_monoid(min, identity));
  point max_corner = reduce(P, make_monoid(max, identity));
  double size = (max_corner-min_corner).Length();
  double stretch = 10.0;
  double radius = stretch*size;
  point center = max_corner + (max_corner-min_corner)/2.0;
  double pi = 3.14159;

  // Generate the bounding points on a circle far outside the bounding box
  for (size_t i=0; i < boundary_size; i++) {
    double x = radius * cos(2*pi*((float) i)/((float) boundary_size));
    double y = radius * sin(2*pi*((float) i)/((float) boundary_size));
    point pt = center + vect(x,y);
    V[i+n] = vertex_t(pt, i + n);
  }

  // Fill with triangles (boundary_size - 2 total)
  simplex_t s = simplex_t(&V[0+n], &V[1+n], &V[2+n], &T[0 + 2*n]); 
  for (size_t i = 3; i < boundary_size; i++)
    s = s.extend(&V[i+n], &T[i - 2 + 2*n]); 
  //return s;
}


// *************************************************************
//    MAIN LOOP
// *************************************************************

void incrementally_add_points(sequence<vertex_t*> v, vertex_t* start) {
  size_t n = v.size();
  
  // various structures needed for each parallel insertion
  size_t max_block_size = (size_t) (n/1000) + 1; // maximum number to try in parallel ??
									  
  sequence<vertex_t*> done(n);  // holds all completed vertices
  sequence<vertex_t*> buffer(max_block_size);// initially empty, holds leftofvers from prev round
  sequence<vertex_t*> remain;  // holds remaining from previous round
  sequence<simplex_t> t(max_block_size);
  sequence<bool> flags(max_block_size);
  auto VQ = tabulate(max_block_size, [&] (size_t i) -> Qs_t {return Qs_t();});
  
  // create a point location structure
  using KNN = k_nearest_neighbors<vertex_t,1>;
  sequence<vertex_t*> init(1,start);
  KNN knn = KNN(init);

  size_t num_done = 0;
  size_t rounds = 0;
  size_t num_failed = 0;
  size_t num_remain = 0;
  size_t num_next_rebuild = 100;
  size_t multiplier = 10;

  while (num_done < n) {
    //if (rounds > 3) abort();

    // every once in a while create a new point location
    // structure using all points inserted so far
    if (num_done >= num_next_rebuild && num_done <= n/multiplier) {
      // cout << "size = " << num_done << endl;
      auto vtxs = parlay::to_sequence(done.cut(0,num_done));
      knn = KNN(vtxs); // should change to pass slice
      num_next_rebuild *= multiplier;
    }

    // determine how many vertices to try in parallel
    size_t num_round = std::min(std::min(1 + num_done/50, n-num_done), max_block_size);
    // 50 is pulled out of a hat
    // cout << "enter loop: " << rounds << ", " << num_round << ", " << num_done << ", " << num_remain << endl;
    
    // for trial vertices find containing triangle, determine cavity 
    // and reserve vertices on boundary of cavity
    parallel_for (0, num_round, [&] (size_t j) {
      buffer[j] = (j < num_remain) ? remain[j] : v[j + num_done];
      vertex_t *u = knn.nearest(buffer[j]);
      t[j] = find(buffer[j], simplex(u->t, 0));
      reserve_for_insert(buffer[j], t[j], &VQ[j]);});
    
    // For trial vertices check if they own their boundary and
    // update mesh if so.  flags[i] is 1 if failed (need to retry)
    parallel_for (0, num_round, [&] (size_t j) {
				  flags[j] = insert(buffer[j], t[j], &VQ[j]);});
    //cout << "here 3: " << flags[0] << endl;

    // Pack failed vertices back onto Q and successful
    // ones up above (needed for point location structure)
    remain = pack(buffer.cut(0,num_round), flags.cut(0,num_round));
    num_remain = remain.size();
    size_t num_done_in_round = num_round - num_remain;
    //cout << "finished " << num_done_in_round << " in round " << rounds << ", " << num_remain << endl;
    auto not_flags = delayed_seq<bool>(num_round, [&] (size_t i) -> bool {return !flags[i];});
    //auto not_flags = tabulate(num_round, [&] (size_t i) -> bool {return !flags[i];});
    pack_out(buffer.cut(0,num_round), not_flags, done.cut(num_done, num_done + num_done_in_round));

    num_failed += num_remain;
    num_done += num_done_in_round;
    rounds++;
  }

  //cout << "n=" << n << "  Total retries=" << failed
  //     << "  Total rounds=" << rounds << endl;
}


// *************************************************************
//    DRIVER
// *************************************************************

triangles<point> delaunay(sequence<point> &P) {
  timer t("delaunay", false);
  t.start();
  size_t boundary_size = 10;
  size_t n = P.size();

  // All vertices needed
  size_t num_vertices = n + boundary_size;
  auto Vertices = sequence<vertex_t>(num_vertices);

  // All triangles needed
  size_t boundary_triangles = (boundary_size - 2);
  size_t num_triangles = 2 * n + boundary_triangles;
  auto Triangles = sequence<triang_t>(num_triangles); 

  // random permutation to put points in a random order
  sequence<size_t> perm = random_permutation<size_t>(n);
  parallel_for(0, n, [&] (size_t i) {
    Vertices[perm[i]] = vertex_t(P[i], i);});

  // give two triangles to each non-boundary vertex
  parallel_for (0, n, [&] (size_t i) {
    Vertices[i].t = &Triangles[2*i];});
  
  // generate boundary points and fill with simplices
  // The boundary points and simplices go at the end,
  // starting at n of Vertices, and 2n of Triangles
  generate_boundary(P, boundary_size, Vertices, Triangles);

  // pointers to first n vertices
  auto V = tabulate(n, [&] (size_t i) -> vertex_t* {
			 return &Vertices[i];});
  vertex_t* v0 = &Vertices[n];
  
  t.next("initialize");
  // main loop to add all points

  incrementally_add_points(V, v0);
  t.next("add points");

  if (CHECK) check_delaunay(Triangles, boundary_size);

  // just the three corner ids for each triangle
  auto result_triangles = tabulate(num_triangles, [&] (size_t i) -> tri {
    vertex_t** vtx = Triangles[i].vtx;
    tri r = {(int) vtx[0]->id, (int) vtx[1]->id, (int) vtx[2]->id};
    return r;});

  // just the points, including the added boundary points
  auto result_points = tabulate(num_vertices, [&] (size_t i) {
    point r = (i < n) ? P[i] : Vertices[i].pt;
    //cout << r[0] << ", " << r[1] << endl;
    return r;});

  t.next("generate output");

  return triangles<point>(result_points, result_triangles);
}
