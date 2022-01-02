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

#include <iostream>
#include <vector>
#include "parlay/primitives.h"
#include "parlay/hash_table.h"
#include "common/get_time.h"
#include "common/topology.h"
#include "common/atomics.h"
#include "refine.h"
#include "common/topology_from_triangles.h"

using std::cout;
using std::endl;
using std::min;
using parlay::hash64;
using parlay::parallel_for;
using parlay::sequence;
using parlay::pack;
using parlay::pack_index;

using vertex_t = vertex<point>;
using simplex_t = simplex<point>;
using triang_t = triangle<point>;
using vect = typename point::vector;

struct Qs {
  vector<vertex<point>*> vertexQ;
  vector<simplex<point>> simplexQ;
  Qs() {
    vertexQ.reserve(50);
    simplexQ.reserve(50);
    // check that actually called
  }
};

using vertexQs = sequence<Qs>;

// *************************************************************
//   PARALLEL HASH TABLE TO STORE WORK QUEUE OF SKINNY TRIANGLES
// *************************************************************

struct hashTriangles {
  typedef triang_t* eType;
  typedef triang_t* kType;
  eType empty() {return NULL;}
  kType getKey(eType v) { return v;}
  size_t hash(kType s) { return hash64(s->id); }
  int cmp(kType s, kType s2) {
    return (s->id > s2->id) ? 1 : ((s->id == s2->id) ? 0 : -1);
  }
  bool cas(eType* p, eType o, eType n) {
    return pbbs::atomic_compare_and_swap(p, o, n);
  }
  bool replaceQ(eType s, eType s2) {return 0;}
};

typedef hashtable<hashTriangles> TriangleTable;
TriangleTable makeTriangleTable(size_t m) {
  return TriangleTable(m,hashTriangles());}

// *************************************************************
//   THESE ARE TAKEN FROM delaunay.C
//   Perhaps should be #included
// *************************************************************

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
void findCavity(simplex_t t, vertex_t *p, Qs *q) {
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
void reserve_for_insert(vertex_t *v, simplex_t t, Qs *q) {
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

// *************************************************************
//   DEALING WITH THE CAVITY
// *************************************************************

inline bool skinnyTriangle(triang_t *t) {
  double minAngle = 30;
  if (minAngleCheck(t->vtx[0]->pt, t->vtx[1]->pt, t->vtx[2]->pt, minAngle))
    return 1;
  return 0;
}

inline bool obtuse(simplex_t t) {
  int o = t.o;
  point p0 = t.t->vtx[(o+1)%3]->pt;
  vect v1 = t.t->vtx[o]->pt - p0;
  vect v2 = t.t->vtx[(o+2)%3]->pt - p0;
  return (v1.dot(v2) < 0.0);
}

inline point circumcenter(simplex_t t) {
  if (t.isTriangle())
    return triangleCircumcenter(t.t->vtx[0]->pt, t.t->vtx[1]->pt, t.t->vtx[2]->pt);
  else { // t.isBoundary()
    point p0 = t.t->vtx[(t.o+2)%3]->pt;
    point p1 = t.t->vtx[t.o]->pt;
    return p0 + (p1-p0)/2.0;
  }
}

// this side affects the simplex_t by moving it into the right orientation
// and setting the boundary if the circumcenter encroaches on a boundary
inline bool checkEncroached(simplex_t& t) {
  if (t.isBoundary()) return 0;
  int i;
  for (i=0; i < 3; i++) {
    if (t.across().isBoundary() && (t.farAngle() > 45.0)) break;
    t = t.rotClockwise();
  }
  if (i < 3) return t.boundary = 1;
  else return 0;
}

bool findAndReserveCavity(vertex_t* v, simplex_t& t, Qs* q) {
  t = simplex_t(v->badT,0);
  if (t.t == NULL) {cout << "refine: nothing in badT" << endl; abort();}
  if (t.t->bad == 0) return 0;

  // if there is an obtuse angle then move across to opposite triangle, repeat
  if (obtuse(t)) t = t.across();
  while (t.isTriangle()) {
    int i;
    for (i=0; i < 2; i++) {
      t = t.rotClockwise();
      if (obtuse(t)) { t = t.across(); break; } 
    }
    if (i==2) break;
  }

  // if encroaching on boundary, move to boundary
  checkEncroached(t);

  // use circumcenter to add (if it is a boundary then its middle)
  v->pt = circumcenter(t);
  reserve_for_insert(v, t, q);
  return 1;
}

// checks if v "won" on all adjacent vertices and inserts point if so
// returns true if "won" and cavity was updated
bool addCavity(vertex_t *v, simplex_t t, Qs *q, TriangleTable& TT) {
  bool flag = 1;
  for (size_t i = 0; i < q->vertexQ.size(); i++) {
    vertex_t* u = (q->vertexQ)[i];
    if (u->reserve == v->id) u->reserve = -1; // reset to -1
    else flag = 0; // someone else with higher priority reserved u
  }
  if (flag) {
    triang_t* t0 = t.t;
    triang_t* t1 = v->t;  // the memory for the two new triang_tangles
    triang_t* t2 = t1 + 1;  
    t1->initialized = 1;
    if (t.isBoundary()) t.splitBoundary(v, t1);
    else {
      t2->initialized = 1;
      t.split(v, t1, t2);
    }

    // update the cavity
    for (size_t i = 0; i<q->simplexQ.size(); i++) 
      (q->simplexQ)[i].flip();
    q->simplexQ.push_back(simplex_t(t0,0));
    q->simplexQ.push_back(simplex_t(t1,0));
    if (!t.isBoundary()) q->simplexQ.push_back(simplex_t(t2,0));

    for (size_t i = 0; i<q->simplexQ.size(); i++) {
      triang_t* t = (q->simplexQ)[i].t;
      if (skinnyTriangle(t)) {
	TT.insert(t); 
	t->bad = 1;}
      else t->bad = 0;
    }
    v->badT = NULL;
  } 
  q->simplexQ.clear();
  q->vertexQ.clear();
  return flag;
}


// *************************************************************
//    MAIN REFINEMENT LOOP
// *************************************************************

// Insert a set of vertices to refine the mesh 
// TT is an initially empty table used to store all the bad
// triangles that are created when inserting vertices
template <typename Slice>
size_t addRefiningVertices(Slice &V, TriangleTable &TT, vertexQs& VQ) {
  size_t n = V.size();
  size_t size = min(VQ.size(), n);
  
  sequence<simplex_t> t(size);
  sequence<bool> flags(size);
  
  size_t top = n;
  size_t num_failed = 0;

  // process all vertices starting just below the top
  while(top > 0) {
    size_t cnt = min<size_t>(size, top);
    size_t offset = top-cnt;

    parallel_for (0, cnt, [&] (size_t j) {
      flags[j] = findAndReserveCavity(V[j+offset], t[j], &VQ[j]);});

    parallel_for (0, cnt, [&] (size_t j) {
      flags[j] = flags[j] && !addCavity(V[j+offset], t[j], &VQ[j], TT);});

    // Pack the failed vertices back onto Q
    auto remain = pack(V.cut(offset,offset+cnt), flags.cut(0,cnt));
    parallel_for (0, remain.size(), [&] (size_t j) {V[j+offset] = remain[j];});
    num_failed += remain.size();
    top = top-cnt+remain.size(); // adjust top, accounting for failed vertices
  }
  return num_failed;
}

// *************************************************************
//    DRIVER
// *************************************************************

#define QSIZE 20000

triangles<point> refineInternal(triangles<point>& Tri) {
  timer t("Delaunay Refine");
  int expandFactor = 4;
  size_t n = Tri.numPoints();
  size_t m = Tri.numTriangles();
  size_t extraVertices = expandFactor*n;
  size_t totalVertices = n + extraVertices;
  size_t totalTriangles = m + 2 * extraVertices;

  sequence<vertex_t> Vertices;
  sequence<triang_t> Triangles;
  sequence<vertex_t*> V(extraVertices);
  
  std::tie(Triangles, Vertices) = topology_from_triangles(Tri, extraVertices);
  t.next("from Triangles");
  
  //  set up extra triangles
  parallel_for (m, totalTriangles, [&] (size_t i) {
    Triangles[i].id = i;
    Triangles[i].initialized = 0;
  });

  //  set up extra vertices
  parallel_for (0, extraVertices, [&] (size_t i) {
    V[i] = new (&Vertices[i+n]) vertex_t(point(0,0), i+n);
    // give each one a pointer to two triangles to use
    V[i]->t = &Triangles[m + 2*i];
  });
  t.next("initializing");

  // these will increase as more are added
  size_t numTriangs = m;
  size_t numPoints = n;

  TriangleTable workQ = makeTriangleTable(numTriangs);
  parallel_for(0, numTriangs, [&] (size_t i) {
    if (skinnyTriangle(&Triangles[i])) {
      workQ.insert(&Triangles[i]);
      Triangles[i].bad = 1;
    }
  });

  vertexQs VQ(QSIZE);
  t.next("Start");

  // Each iteration processes all bad triangles from the workQ while
  // adding new bad triangles to a new queue
  while (1) {
    sequence<triang_t*> badTT = workQ.entries();

    // packs out triangles that are no longer bad
    auto flags = tabulate(badTT.size(), [&] (size_t i) -> bool {
      return badTT[i]->bad;});
    auto badT = pack(badTT, flags);
    size_t numBad = badT.size();

    cout << "numBad = " << numBad << endl;
    if (numBad == 0) break;
    if (numPoints + numBad > totalVertices) {
      cout << "ran out of vertices" << endl;
      abort();
    }
    size_t offset = numPoints - n;

    // allocate 1 vertex per bad triangle and assign triangle to it
    parallel_for (0, numBad, [&] (size_t i) {
      badT[i]->bad = 2; // used to detect whether touched
      V[i + offset]->badT = badT[i];
    });

    // the new empty work queue
    workQ = makeTriangleTable(numBad);

    // This does all the work adding new vertices, and any new bad
    // triangles to the workQ
    auto Vtx = V.cut(offset, offset+numBad);
    addRefiningVertices(Vtx, workQ, VQ);

    // push any bad triangles that were left untouched onto the Q
    parallel_for (0, numBad, [&] (size_t i) {
      if (badT[i]->bad==2) workQ.insert(badT[i]);});

    numPoints += numBad;
    numTriangs += 2*numBad;
  }

  t.next("refinement");
  std::cout << numTriangs << " : " << Vertices.size() << " : " << numPoints << std::endl;
  
  // Extract Vertices for result
  auto flag = tabulate(numPoints, [&] (size_t i) -> bool {
    return (Vertices[i].badT == NULL);});

  std::cout << "here" << std::endl;
  sequence<size_t> I = pack_index(flag);
  size_t n0 = I.size();
  sequence<point> rp(n0);

  std::cout << "here2" << std::endl;
  parallel_for (0, n0, [&] (size_t i) {
    Vertices[I[i]].id = i;
    rp[i] = Vertices[I[i]].pt;
  });
  cout << "total points = " << n0 << endl;

  // Extract Triangles for result
  I = pack_index(tabulate(numTriangs, [&] (size_t i) -> bool {
	 return Triangles[i].initialized;}));
							  
  auto rt = tabulate(I.size(), [&] (size_t i) -> tri {
    auto t = Triangles[I[i]];
    tri r = {t.vtx[0]->id, t.vtx[1]->id, t.vtx[2]->id};
    return r;});

  cout << "total triangles = " << I.size() << endl;
  t.next("finish");
  return triangles<point>(std::move(rp), std::move(rt));
}

triangles<point> refine(triangles<point> &Tri) {
  return refineInternal(Tri);
}
