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
// "Multipole" and "Local" expansion and translations between them.
// For the expansions it uses a modified version of the multipole
// translation code from the PETFMM library using spherical
// harmonics.  The translations are implemented in spherical.h and can
// be changed for any other routines that support the public interface
// of the transform structure.
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
//   BOXSIZE -- The max number of particles in each leaf of the tree

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

#define CHECK 0

// Following for 1e-1 accuracy (1.05 seconds for 1million 8 cores)
//#define ALPHA 1.9
//#define terms 3
//#define BOXSIZE 25

// Following for 1e-3 accuracy (4 seconds for 1million 8 cores)
//#define ALPHA 2.2
//#define terms 7
//#define BOXSIZE 60

// Following for 1e-6 accuracy (12.5 seconds for 1million insphere 8 cores)
#define ALPHA 2.65
#define terms 12
#define BOXSIZE 250 // 130

// Following for 1e-9 accuracy (40 seconds for 1million 8 cores)
//#define ALPHA 3.0
//#define terms 17
//#define BOXSIZE 250


double check(sequence<particle*> const &p) {
  size_t n = p.size();
  size_t nCheck = min<size_t>(n,200);
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

Transform<terms>* TRglobal = new Transform<terms>();

using box = pair<point,point>;
using vect3d = typename point::vector;

struct node {
  using edge = pair<node*, size_t>;
  node* left;
  node* right;
  sequence<particle*> particles;
  size_t n;
  box b;
  innerExpansion* InExp;
  outerExpansion* OutExp;
  vector<node*> indirectNeighbors;
  vector<edge> leftNeighbors;
  vector<edge> rightNeighbors;
  sequence<vect3d*> hold;
  bool leaf() {return left == NULL;}
  node() {}
  point center() { return b.first + (b.second-b.first)/2.0;}
  double radius() { return (b.second - b.first).Length()/2.0;}
  double lmax() {
    vect3d d = b.second-b.first;
    return max(d.x,max(d.y,d.z));
  }
  void allocateExpansions() {
    InExp = (innerExpansion*) parlay::p_malloc(sizeof(innerExpansion));
    OutExp = (outerExpansion*) parlay::p_malloc(sizeof(outerExpansion));
    new (InExp) innerExpansion(TRglobal, center());
    new (OutExp) outerExpansion(TRglobal, center());
  }
  node(node* L, node* R, size_t n, box b)
    : left(L), right(R), n(n), b(b) {
    allocateExpansions();
  }
  node(parlay::sequence<particle*> P, box b) 
    : left(NULL), right(NULL), particles(std::move(P)), b(b) {
    n = particles.size();
    allocateExpansions();
  }
};


using edge = pair<node*, size_t>;

template <typename ParticleSlice>
node* buildTree(ParticleSlice particles, ParticleSlice Tmp, size_t effective_size) {
  
  //if (depth > 100) abort();
  size_t n = particles.size();
  size_t en = std::max(effective_size, n);

  auto minmax = [] (box a, box b) {
    return box((a.first).minCoords(b.first),
	       (a.second).maxCoords(b.second));};
  auto pairs = parlay::delayed_map(particles, [&] (particle* p) {
      return box(p->pt, p->pt);});
  box b = parlay::reduce(pairs, parlay::make_monoid(minmax,pairs[0]));
										      
  if (en < BOXSIZE || n < 10) { 
    node* r = (node*) parlay::p_malloc(sizeof(node));
    new (r) node(parlay::to_sequence(particles), b);
    return r;
  }

  size_t d = 0;
  double Delta = 0.0;
  for (int i=0; i < 3; i++) {
    if (b.second[i] - b.first[i] > Delta) {
      d = i;
      Delta = b.second[i] - b.first[i];
    }
  }
  
  double splitpoint = (b.first[d] + b.second[d])/2.0;
  
  auto flagsLeft = parlay::map(particles, [&] (particle* p) -> bool {
      return p->pt[d] < splitpoint;});
  auto flagsRight = parlay::delayed_map(flagsLeft, [] (bool x) {
      return !x;});

  size_t nl = parlay::pack_into_uninitialized(particles, flagsLeft, Tmp);
  parlay::pack_into_uninitialized(particles, flagsRight, Tmp.cut(nl, n));
  parlay::copy(Tmp, particles);
  size_t en_child = .44 * en;

  node* L;
  node* R;
  parlay::par_do(
    [&] () {L = buildTree(particles.cut(0,nl), Tmp.cut(0,nl), en_child);},
    [&] () {R = buildTree(particles.cut(nl,n), Tmp.cut(nl,n), en_child);});

  return new node(L, R, n, b);
}

bool far(node* a, node* b) {
  double rmax = max(a->radius(), b->radius());
  double r = (a->center() - b->center()).Length();
  return r >= (ALPHA * rmax);
}

// used to count the number of interactions
struct ipair {
  long direct;
  long indirect;
  ipair() {}
  ipair(long a, long b) : direct(a), indirect(b) {}
  ipair operator+ (ipair b) {
    return ipair(direct + b.direct, indirect + b.indirect);}
};


ipair interactions(node* Left, node* Right) {
  if (far(Left,Right)) {
    Left->indirectNeighbors.push_back(Right); 
    Right->indirectNeighbors.push_back(Left); 
    return ipair(0,2);
  } else {
    if (!Left->leaf() && (Left->lmax() >= Right->lmax() || Right->leaf())) {
      ipair x = interactions(Left->left, Right);
      ipair y = interactions(Left->right, Right);
      return x + y;
    } else if (!Right->leaf()) {
      ipair x = interactions(Left, Right->left);
      ipair y = interactions(Left, Right->right);
      return x + y;
    } else { // both are leaves
      if (Right->n > Left->n) swap(Right,Left);
      size_t rn = Right->leftNeighbors.size();
      size_t ln = Left->rightNeighbors.size();
      Right->leftNeighbors.push_back(edge(Left,ln)); 
      Left->rightNeighbors.push_back(edge(Right,rn));
      return ipair(Right->n*Left->n,0);
    }
  }
}

// could be parallelized, but so fast it does not help
ipair interactions(node* tr) {
  if (!tr->leaf()) {
    ipair x, y, z; 
    x = interactions(tr->left);
    y = interactions(tr->right);
    z = interactions(tr->left,tr->right);
    return x + y + z;
  } else return ipair(0,0);
}

size_t numLeaves(node* tr) {
  if (tr->leaf()) return 1;
  else return(numLeaves(tr->left)+numLeaves(tr->right));
}

// Translate from multipole to local expansion along all far-field
// interactions.
void doIndirect(node* tr) {
  //tr->OutExp = new outerExpansion(TRglobal, tr->center());
  for (size_t i = 0; i < tr->indirectNeighbors.size(); i++) 
    tr->OutExp->addTo(tr->indirectNeighbors[i]->InExp);
  if (!tr->leaf()) {
    parlay::par_do([&] () {doIndirect(tr->left);},
		   [&] () {doIndirect(tr->right);});
  }
}

// Translate and accumulate multipole expansions up the tree,
// including translating particles to expansions at the leaves.
void upSweep(node* tr) {
  //tr->InExp = new innerExpansion(TRglobal, tr->center());
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

// Translate and accumulate local expansions down the tree,
// including applying them to all particles at the leaves.
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

// Calculates the direct forces between all pairs of particles in two nodes.
// Directly updates forces in Left, and places forces for ngh in hold
// This avoid a race condition on modifying ngh while someone else is
// updating it.
void direct(node* Left, node* ngh, vect3d* hold) {
  particle** LP = (Left->particles).data();
  particle** RP = (ngh->particles).data();
  size_t nl = Left->n;
  size_t nr = ngh->n;
  vect3d holdLeft[nr];
  for (size_t j=0; j < nr; j++) 
    holdLeft[j] = vect3d(0.,0.,0.);
  for (size_t i=0; i < nl; i++) {
    vect3d frc(0.,0.,0.);
    particle* pa = LP[i];
    for (size_t j=0; j < nr; j++) {
      particle* pb = RP[j];
      vect3d v = (pb->pt) - (pa->pt);
      double r2 = v.dot(v);
      vect3d force = (v * (pa->mass * pb->mass / (r2*sqrt(r2))));;
      holdLeft[j] = holdLeft[j] - force;
      frc = frc + force;
    }
    pa->force = pa->force + frc;
  }
  for (size_t j=0; j < nr; j++) 
    hold[j] = holdLeft[j];
}

void self(node* Tr) {
  particle** PP = (Tr->particles).data();
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

// Calculates the direct interactions between and within leaves.
// Since the forces are symmetric, this calculates the force on one
// side (rightNeighbors) while storing them away (in hold).
// It then goes over the other side (leftNeighbors) picking up
// the precalculated results (from hold).
// It does not update both sides immediately since that would 
// generate a race condition.
void doDirect(node* a) {
  size_t nleaves = numLeaves(a);
  sequence <node*> Leaves(nleaves);
  getLeaves(a, Leaves.data());

  sequence<size_t> counts(nleaves+1);
  
  // For node i in Leaves, counts[i] will contain the total number 
  // of its direct interactions.
  parlay::parallel_for (0, nleaves, [&] (size_t i) {
    counts[i] = 0;
    for (size_t j = 0; j < Leaves[i]->rightNeighbors.size(); j++)
      counts[i] += Leaves[i]->rightNeighbors[j].first->n;
    });

  // The following allocates space for "hold" avoiding a malloc.
  counts[nleaves] = 0;
  size_t total = parlay::scan_inplace(counts);
  sequence<vect3d> hold(total);
  
  // calculates interactions and neighbors results in hold
  parlay::parallel_for (0, nleaves, [&] (size_t i) {
    vect3d* lhold = hold.begin() + counts[i];
    size_t num_ngh = Leaves[i]->rightNeighbors.size();
    Leaves[i]->hold = sequence<vect3d*>(num_ngh);
    for (size_t j =0; j < num_ngh; j++) {
      Leaves[i]->hold[j] = lhold;
      node* ngh = Leaves[i]->rightNeighbors[j].first;
      direct(Leaves[i], ngh, lhold);
      lhold += ngh->n;
    }
  }, 1); // 1 indicates granularity of 1
  
  // picks up results from neighbors that were left in hold
  parlay::parallel_for (0, nleaves, [&] (size_t i) {
    for (size_t j = 0; j < Leaves[i]->leftNeighbors.size(); j++) {
      node* L = Leaves[i];
      edge e = L->leftNeighbors[j];
      vect3d* hold = e.first->hold[e.second];
      for (size_t k=0; k < Leaves[i]->n; k++) 
	L->particles[k]->force = L->particles[k]->force + hold[k];
    }
    });

  // calculate forces within a node
  parlay::parallel_for (0, nleaves, [&] (size_t i) {self(Leaves[i]);});
}

// *************************************************************
//   STEP
// *************************************************************

// takes one step and places forces in particles[i]->force
void stepBH(sequence<particle*> &particles, double alpha) {
  timer t("CK nbody");
  size_t n = particles.size();
  TRglobal->precompute();

  parlay::parallel_for (0, n, [&] (size_t i) {
      particles[i]->force = vect3d(0.,0.,0.);});

  sequence<particle*> Tmp(n);
  sequence<particle*> part_copy = particles;

  // build the CK tree
  node* a = buildTree(parlay::make_slice(part_copy), parlay::make_slice(Tmp), 0);
  t.next("build tree");

  // Sweep up the tree calculating multipole expansions for each node
  upSweep(a);
  t.next("up sweep");

  // Determine all far-field interactions using the CK method
  ipair z = interactions(a);
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
