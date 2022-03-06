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

// This is an implementation of the Callahan-Kosaraju (CK) algorithm
// for n-body simulation.
//
//   Paul Callahan and S. Rao Kosaraju
//   A decomposition of multi-dimensional point-sets with applications 
//   to k-nearest-neighbors and n-body potential fields 
//   ACM Symposium on Theory of Computation, 1992
// 
// It uses similar ideas to the Greengard-Rothkin FMM method but is
// more flexible for umbalanced trees.  As with FMM it uses
// "Multipole" (or "inner") and "Local" (or "outer") expansion and
// translations between them.  For the expansions it uses a modified
// version of the multipole translation code from the PETFMM library
// using spherical harmonics.  The translations are implemented in
// spherical.h and can be changed for any other routines that support
// the public interface of the transform structure.
//
// Similarly to most FMM-based codes it works in the following steps
//   1) build the CK tree recursively (similar to a k-d tree)
//   2) calculate multipole expansions going up the tree
//   3) figure out all far-field interactions using the CK method
//   4) translate all multipole to local expansions along 
//      the far-field interactions calculated in 3).
//   5) propagate local expansions down the tree
//   6) finally add in all direct leaf-leaf interactions
//
// The accuracy can be adjusted using the following parameters
//   ALPHA -- Controls distance which is considered far-field.
//            It is the min ratio of the larger of two interacting 
//            boxes to the distance between them
//   terms -- Number of terms in the expansions
// The performance can be adjusted with
//   BOXSIZE -- The max number of particles in each leaf of the tree
//      this also slightly affects accuracy (smaller is better)

#include <iostream>
#include <vector>
#include "common/get_time.h"
#include "common/geometry.h"
#include "parlay/primitives.h"
#include "spherical.h"
#include "nbody.h"

using namespace std;
using parlay::sequence;
using vect3d = vect;

#define CHECK 1

// Following for 1e-3 accuracy
//#define ALPHA 2.2
//#define terms 7
//#define BOXSIZE 150

// Following for 1e-6 accuracy (2.5x slower than above)
#define ALPHA 2.6
#define terms 12  
#define BOXSIZE 250

// Following for 1e-9 accuracy (2.2x slower than above)
// #define ALPHA 3.0
// #define terms 17
// #define BOXSIZE 550

// Following for 1e-12 accuracy (1.8x slower than above)
//#define ALPHA 3.2
//#define terms 22
//#define BOXSIZE 700

double check(sequence<particle*> const &p) {
  size_t n = p.size();
  size_t nCheck = min<size_t>(n, 200);
  sequence<double> Err(nCheck);
  
  parlay::parallel_for (0, nCheck, [&] (size_t i) {
    size_t idx = parlay::hash64(i)%n;
    vect3d force(0.,0.,0.);
    for (size_t j=0; j < n; j++) {
      if (idx != j) {
	vect3d v = (p[j]->pt) - (p[idx]->pt);
	double r2 = v.dot(v);
	force = force + (v * (p[j]->mass * p[idx]->mass / (r2*sqrt(r2))));
      }
    }
    Err[i] = (force - p[idx]->force).Length()/force.Length();
    });
  double total = 0.0;
  for(int i=0; i < nCheck; i++) 
    total += Err[i];
  return total/nCheck;
}

// *************************************************************
//    FORCE CALCULATIONS
// *************************************************************

// *************************************************************
//  Inner expansions (also called multipole expansion)
//  The spherical harmonic expansion of a set of nearby points around
//  a center for estimating forces at a distance.
// *************************************************************
struct innerExpansion {
  Transform<terms>* TR;
  complex<double> coefficients[terms*terms];
  point center;
  void addTo(point pt, double mass) {
    TR->P2Madd(coefficients, mass, center, pt);
  }
  void addTo(innerExpansion* y) {
    TR->M2Madd(coefficients, center, y->coefficients, y->center);
  }
  innerExpansion(Transform<terms>* _TR, point _center) : TR(_TR), center(_center) {
    for (size_t i=0; i < terms*terms; i++) coefficients[i] = 0.0;
  }
  vect3d force(point y, double mass) {
    vect3d result;
    double potential;
    TR->M2P(potential, result, y, coefficients, center);
    result = result*mass;
    return result;
  }
  innerExpansion() {}
};

parlay::type_allocator<innerExpansion> inner_pool;

// *************************************************************
//  Outer expansions (also called local)
//  The inverse spherical harmonic expansion of a set of distant
//  points around a center for estimating forces for nearby points.
// *************************************************************
struct outerExpansion {
  Transform<terms>* TR;
  complex<double> coefficients[terms*terms];
  point center;
  void addTo(innerExpansion* y) {
    TR->M2Ladd(coefficients, center, y->coefficients, y->center);}
  void addTo(outerExpansion* y) {
    TR->L2Ladd(coefficients, center, y->coefficients, y->center);
  }
  vect3d force(point y, double mass) {
    vect3d result;
    double potential;
    TR->L2P(potential, result, y, coefficients, center);
    result = result*mass;
    return result;
  }
  outerExpansion(Transform<terms>* _TR, point _center) : TR(_TR), center(_center) {
    for (size_t i=0; i < terms*terms; i++) coefficients[i] = 0.0;
  }
  outerExpansion() {}
};

parlay::type_allocator<outerExpansion> outer_pool;

// Set global constants for spherical harmonics
Transform<terms>* TRglobal = new Transform<terms>();

using box = pair<point,point>;
using vect3d = typename point::vector;

// *************************************************************
//  A node in the CK tree
//  Either a leaf (if children are null) or internal node.
//  If a leaf contains a set of points
//  If an internal node contains a left and right child as is
//  augmented with first inner than outer expansions.
//  The leftNeighbors and rightNeighbors contain edges in the CK
//  well separated decomposition.
// *************************************************************
struct node {
  using edge = pair<node*, size_t>;
  node* left;
  node* right;
  sequence<particle*> particles;
  sequence<particle> particles_d;
  size_t n;
  box b;
  innerExpansion* InExp;
  outerExpansion* OutExp;
  vector<node*> indirectNeighbors;
  vector<edge> leftNeighbors;
  vector<edge> rightNeighbors;
  sequence<sequence<vect3d>> hold;
  bool leaf() {return left == NULL;}
  node() {}
  point center() { return b.first + (b.second-b.first)/2.0;}
  double radius() { return (b.second - b.first).Length()/2.0;}
  double lmax() {
    vect3d d = b.second-b.first;
    return max(d.x,max(d.y,d.z));
  }
  void allocateExpansions() {
    InExp = inner_pool.allocate(TRglobal, center());
    OutExp = outer_pool.allocate(TRglobal, center());
  }
  node(node* L, node* R, size_t n, box b)
    : left(L), right(R), n(n), b(b) {
    allocateExpansions();
  }
  node(parlay::sequence<particle*> P, box b) 
    : left(NULL), right(NULL), particles(std::move(P)), b(b) {
    n = particles.size();
    particles_d = parlay::map(particles, [] (auto p) {return *p;});
    allocateExpansions();
  }
};

size_t numLeaves(node* tr) {
  if (tr->leaf()) return 1;
  else return(numLeaves(tr->left)+numLeaves(tr->right));
}

parlay::type_allocator<node> node_pool;

using edge = pair<node*, size_t>;

// *************************************************************
//  Build the CK tree
//  Similar to a kd-tree but always split along widest dimension
//  of the points instead of the next round-robin dimension.
// *************************************************************
template <typename Particles>
node* buildTree(Particles& particles, size_t effective_size) {
  
  size_t n = particles.size();
  size_t en = std::max(effective_size, n);

  auto minmax = [] (box a, box b) {
    return box((a.first).minCoords(b.first),
	       (a.second).maxCoords(b.second));};
  auto pairs = parlay::delayed_map(particles, [&] (particle* p) {
      return box(p->pt, p->pt);});
  box b = parlay::reduce(pairs, parlay::make_monoid(minmax,pairs[0]));
										      
  if (en < BOXSIZE || n < 10) 
    return node_pool.allocate(parlay::to_sequence(particles), b);

  size_t d = 0;
  double Delta = 0.0;
  for (int i=0; i < 3; i++) {
    if (b.second[i] - b.first[i] > Delta) {
      d = i;
      Delta = b.second[i] - b.first[i];
    }
  }
  
  double splitpoint = (b.first[d] + b.second[d])/2.0;

  auto isLeft = parlay::delayed_map(particles, [&] (particle* p) {
      return std::pair(p->pt[d] < splitpoint, p);});
  auto foo = parlay::group_by_index(isLeft, 2);
  particles.clear();

  auto r = parlay::map(foo, [&] (auto& x) {
      return buildTree(x, .4 * en);}, 1);
  return node_pool.allocate(r[0], r[1], n, b);
}


// *************************************************************
//  Determine if a point is far enough to use approximation.
// *************************************************************
bool far(node* a, node* b) {
  double rmax = max(a->radius(), b->radius());
  double r = (a->center() - b->center()).Length();
  return r >= (ALPHA * rmax);
}

// *************************************************************
// Used to count the number of interactions, just for performance
// statistics not needed for correctness.
// *************************************************************
struct interactions_count {
  long direct;
  long indirect;
  interactions_count() {}
  interactions_count(long a, long b) : direct(a), indirect(b) {}
  interactions_count operator+ (interactions_count b) {
    return interactions_count(direct + b.direct, indirect + b.indirect);}
};

// *************************************************************
// The following two functions are the core of the CK method.
// They calculate the "well separated decomposition" of the points.
// *************************************************************
interactions_count interactions(node* Left, node* Right) {
  if (far(Left,Right)) {
    Left->indirectNeighbors.push_back(Right); 
    Right->indirectNeighbors.push_back(Left); 
    return interactions_count(0,2);
  } else {
    if (!Left->leaf() && (Left->lmax() >= Right->lmax() || Right->leaf())) {
      interactions_count x = interactions(Left->left, Right);
      interactions_count y = interactions(Left->right, Right);
      return x + y;
    } else if (!Right->leaf()) {
      interactions_count x = interactions(Left, Right->left);
      interactions_count y = interactions(Left, Right->right);
      return x + y;
    } else { // both are leaves
      if (Right->n > Left->n) swap(Right,Left);
      size_t rn = Right->leftNeighbors.size();
      size_t ln = Left->rightNeighbors.size();
      Right->leftNeighbors.push_back(edge(Left,ln)); 
      Left->rightNeighbors.push_back(edge(Right,rn));
      return interactions_count(Right->n*Left->n,0);
    }
  }
}

// Could be parallelized but would require avoiding push_back.
// Currently not a bottleneck so left serial.
interactions_count interactions(node* tr) {
  if (!tr->leaf()) {
    interactions_count x, y, z; 
    x = interactions(tr->left);
    y = interactions(tr->right);
    z = interactions(tr->left,tr->right);
    return x + y + z;
  } else return interactions_count(0,0);
}

// *************************************************************
// Translate from inner (multipole) expansion to outer (local)
// expansion along all far-field interactions.
// *************************************************************
void doIndirect(node* tr) {
  for (size_t i = 0; i < tr->indirectNeighbors.size(); i++) 
    tr->OutExp->addTo(tr->indirectNeighbors[i]->InExp);
  if (!tr->leaf()) {
    parlay::par_do([&] () {doIndirect(tr->left);},
		   [&] () {doIndirect(tr->right);});
  }
}

// *************************************************************
// Translate and accumulate inner (multipole) expansions up the tree,
// including translating particles to expansions at the leaves.
// *************************************************************
void upSweep(node* tr) {
  if (tr->leaf()) {
    for (size_t i=0; i < tr->n; i++) {
      particle* P = tr->particles[i];
      tr->InExp->addTo(P->pt, P->mass);
    }
  } else {
    parlay::par_do([&] () {upSweep(tr->left);},
		   [&] () {upSweep(tr->right);});
    tr->InExp->addTo(tr->left->InExp);
    tr->InExp->addTo(tr->right->InExp);
  }
}

// *************************************************************
// Translate and accumulate outer (local) expansions down the tree,
// including applying them to all particles at the leaves.
// *************************************************************
void downSweep(node* tr) {
  if (tr->leaf()) {
    for (size_t i=0; i < tr->n; i++) {
      particle* P = tr->particles[i];
      P->force = P->force + tr->OutExp->force(P->pt, P->mass);
    }
  } else {
    parlay::par_do([&] () {tr->left->OutExp->addTo(tr->OutExp);
	                   downSweep(tr->left);},
		   [&] () {tr->right->OutExp->addTo(tr->OutExp);
		           downSweep(tr->right);});
  }
}

// puts the leaves of tree tr into the array Leaves in left to right order
size_t getLeaves(node* tr, node** Leaves) {
  if (tr->leaf()) {
    Leaves[0] = tr;
    return 1;
  } else {
    size_t l = getLeaves(tr->left, Leaves);
    size_t r = getLeaves(tr->right, Leaves + l);
    return l + r;
  }
}

// *************************************************************
// Calculates the direct forces between all pairs of particles in two nodes.
// Directly updates forces in Left, and places forces for ngh in hold
// This avoid a race condition on modifying ngh while someone else is
// updating it.
// *************************************************************
auto direct(node* Left, node* ngh) {
  auto LP = (Left->particles).data();
  auto R = (ngh->particles_d).data();
  size_t nl = Left->n;
  size_t nr = ngh->n;
  parlay::sequence<vect3d> hold(nr, vect3d(0.,0.,0.));
  for (size_t i=0; i < nl; i++) {
    vect3d frc(0.,0.,0.);
    particle pa = *LP[i];
    for (size_t j=0; j < nr; j++) {
      particle& pb = R[j];
      vect3d v = pb.pt - pa.pt;
      double r2 = v.dot(v);
      vect3d force = (v * (pa.mass * pb.mass / (r2*sqrt(r2))));;
      hold[j] = hold[j] - force;
      frc = frc + force;
    }
    LP[i]->force = LP[i]->force + frc;
  }
  return hold;
}

// *************************************************************
// Calculates local forces within a leaf
// *************************************************************
void self(node* Tr) {
  auto PP = (Tr->particles).data();
  for (size_t i=0; i < Tr->n; i++) {
    particle* pa = PP[i];
    for (size_t j=i+1; j < Tr->n; j++) {
	particle* pb = PP[j];
	vect3d v = (pb->pt) - (pa->pt);
	double r2 = v.dot(v);
	vect3d force = (v * (pa->mass * pb->mass / (r2*sqrt(r2))));
	pb->force = pb->force - force;
	pa->force = pa->force + force;
      }
  }
}

// *************************************************************
// Calculates the direct interactions between and within leaves.
// Since the forces are symmetric, this calculates the force on one
// side (rightNeighbors) while storing them away (in hold).
// It then goes over the other side (leftNeighbors) picking up
// the precalculated results (from hold).
// It does not update both sides immediately since that would 
// generate a race condition.
// *************************************************************
void doDirect(node* a) {
  size_t nleaves = numLeaves(a);
  sequence <node*> Leaves(nleaves);
  getLeaves(a, Leaves.data());

  // calculates interactions and put neighbor's results in hold
  parlay::parallel_for (0, nleaves, [&] (size_t i) {
      size_t rn = Leaves[i]->rightNeighbors.size();
      Leaves[i]->hold = parlay::tabulate(rn, [&] (size_t j) {
	  return direct(Leaves[i], Leaves[i]->rightNeighbors[j].first);}, rn);}, 1);

  // picks up results from neighbors that were left in hold
  parlay::parallel_for (0, nleaves, [&] (size_t i) {
    for (size_t j = 0; j < Leaves[i]->leftNeighbors.size(); j++) {
      node* L = Leaves[i];
      auto [u, v] = L->leftNeighbors[j];
      for (size_t k=0; k < Leaves[i]->n; k++) 
	L->particles[k]->force = L->particles[k]->force + u->hold[v][k];
    }}, 1);

  // calculate forces within a node
  parlay::parallel_for (0, nleaves, [&] (size_t i) {self(Leaves[i]);});
}

// *************************************************************
// STEP
// takes one step and places forces in particles[i]->force
// *************************************************************
void stepBH(sequence<particle*> &particles, double alpha) {
  timer t("CK nbody");
  size_t n = particles.size();
  TRglobal->precompute();

  parlay::parallel_for (0, n, [&] (size_t i) {
      particles[i]->force = vect3d(0.,0.,0.);});

  sequence<particle*> part_copy = particles;

  // build the CK tree
  node* a = buildTree(part_copy, 0);
  t.next("build tree");

  // Sweep up the tree calculating multipole expansions for each node
  upSweep(a);
  t.next("up sweep");

  // Determine all far-field interactions using the CK method
  interactions_count z = interactions(a);
  t.next("interactions");

  // Translate multipole to local expansions along the far-field
  // interactions
  doIndirect(a);
  t.next("do Indirect");

  // Translate the local expansions down the tree to the leaves
  downSweep(a);
  t.next("down sweep");

  // Add in all the direct (near-field) interactions
  doDirect(a);
  t.next("do Direct");

  cout << "Direct = " << (long) z.direct << " Indirect = " << z.indirect
       << " Boxes = " << numLeaves(a) << endl;
  if (CHECK) {
    cout << "  Sampled RMS Error = "<< check(particles) << endl;
    t.next("check");
  }
}

void nbody(sequence<particle*> &particles) { 
  stepBH(particles, ALPHA); }
