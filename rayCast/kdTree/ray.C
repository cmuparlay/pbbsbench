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

#include <limits>
#include <algorithm>
#include "common/geometry.h"
#include "common/get_time.h"
#include "parlay/primitives.h"
#include "ray.h"
#include "kdTree.h"
#include "rayTriangleIntersect.h"
using namespace std;

using parlay::tabulate;
using parlay::parallel_for;
using parlay::scan_inplace;
using parlay::sequence;
using parlay::par_do;
using parlay::pack;
using parlay::sort_inplace;
using parlay::to_sequence;

int CHECK = 0;  // if set checks 10 rays against brute force method
int STATS = 0;  // if set prints out some tree statistics

// Constants for deciding when to stop recursion in building the KDTree
float CT = 6.0;
float CL = 1.25;
float maxExpand = 1.6;
int maxRecursionDepth = 25;

// Constant for switching to sequential versions
int minParallelSize = 500000;

using vect = typename point::vector;
using triangs = triangles<point>;
using ray_t = ray<point>;

float boxSurfaceArea(BoundingBox B) {
  float r0 = B[0].max-B[0].min;
  float r1 = B[1].max-B[1].min;
  float r2 = B[2].max-B[2].min;
  return 2*(r0*r1 + r1*r2 + r0*r2);
}

float epsilon = .0000001;
range fixRange(float minv, float maxv) {
  if (minv == maxv) return range(minv,minv+epsilon);
  else return range(minv,maxv);
}

inline float inBox(point p, BoundingBox B) {
  return (p.x >= (B[0].min - epsilon) && p.x <= (B[0].max + epsilon) &&
	  p.y >= (B[1].min - epsilon) && p.y <= (B[1].max + epsilon) &&
	  p.z >= (B[2].min - epsilon) && p.z <= (B[2].max + epsilon));
}

// sequential version of best cut
cutInfo bestCutSerial(sequence<event> const &E, range r, range r1, range r2) {
  double flt_max = std::numeric_limits<double>::max();
  index_t n = E.size();
  if (r.max - r.min == 0.0) return cutInfo(flt_max, r.min, n, n);
  float area = 2 * (r1.max-r1.min) * (r2.max-r2.min);
  float diameter = 2 * ((r1.max-r1.min) + (r2.max-r2.min));

  // calculate cost of each possible split
  index_t inLeft = 0;
  index_t inRight = n/2;
  float minCost = flt_max;
  index_t k = 0;
  index_t rn = inLeft;
  index_t ln = inRight;
  for (index_t i=0; i <n; i++) {
    float cost;
    if (IS_END(E[i])) inRight--;
    float leftLength = E[i].v - r.min;
    float leftArea = area + diameter * leftLength;
    float rightLength = r.max - E[i].v;
    float rightArea = area + diameter * rightLength;
    cost = (leftArea * inLeft + rightArea * inRight);
    if (cost < minCost) {
      rn = inRight;
      ln = inLeft;
      minCost = cost;
      k = i;
    }
    if (IS_START(E[i])) inLeft++;
  }
  return cutInfo(minCost, E[k].v, ln, rn);
}

// parallel version of best cut
cutInfo bestCut(sequence<event> const &E, range r, range r1, range r2) {
  //cout << "bestcut : " << E[0].v << ", " << IS_START(E[0]) << ", " << GET_INDEX(E[0]) << endl;
  index_t n = E.size();
  if (n < minParallelSize)
    return bestCutSerial(E, r, r1, r2);
  double flt_max = std::numeric_limits<double>::max();
  if (r.max - r.min == 0.0) return cutInfo(flt_max, r.min, n, n);

  // area of two orthogonal faces
  float orthogArea = 2 * ((r1.max-r1.min) * (r2.max-r2.min));

  // length of diameter of orthogonal face
  float diameter = 2 * ((r1.max-r1.min) + (r2.max-r2.min));

  // count number that end before i
  auto upperC = tabulate(n, [&] (index_t i) -> index_t {return IS_END(E[i]);});
  scan_inplace(upperC);

  // calculate cost of each possible split location
  auto cost = tabulate(n, [&] (index_t i) -> float {
    index_t inLeft = i - upperC[i];
    index_t inRight = n/2 - (upperC[i] + IS_END(E[i]));
    float leftLength = E[i].v - r.min;
    float leftArea = orthogArea + diameter * leftLength;
    float rightLength = r.max - E[i].v;
    float rightArea = orthogArea + diameter * rightLength;
    return (leftArea * inLeft + rightArea * inRight);
			  });

  // find minimum across all (maxIndex with less is minimum)
  index_t k = parlay::min_element(cost) - cost.begin();

  //cout << "c : " << cost[k] << ", " << orthogArea << ", " << diameter <<  endl;
  float c = cost[k];
  index_t ln = k - upperC[k];
  index_t rn = n/2 - (upperC[k] + IS_END(E[k]));
  return cutInfo(c, E[k].v, ln, rn);
}

using eventsPair = pair<sequence<event>, sequence<event>>;

eventsPair splitEventsSerial(sequence<range> const &boxes,
			     sequence<event> const & events, 
			     float cutOff) {
  index_t l = 0;
  index_t r = 0;
  index_t n = events.size();
  sequence<event> eventsLeft(n);
  sequence<event> eventsRight(n);
  for (index_t i=0; i < n; i++) {
    index_t b = GET_INDEX(events[i]);
    if (boxes[b].min < cutOff) {
      eventsLeft[l++] = events[i];
      if (boxes[b].max > cutOff) 
	eventsRight[r++] = events[i]; 
    } else eventsRight[r++] = events[i]; 
  }
  return eventsPair(to_sequence(eventsLeft.cut(0,l)),
		    to_sequence(eventsRight.cut(0,r)));
}

eventsPair splitEvents(sequence<range> const &boxes,
		       sequence<event> const &events, 
		       float cutOff) {
  index_t n = events.size();
  if (n < minParallelSize)
    return splitEventsSerial(boxes, events, cutOff);
  auto lower = sequence<bool>::uninitialized(n);
  auto upper = sequence<bool>::uninitialized(n);

  parallel_for (0, n, [&] (index_t i) {
    index_t b = GET_INDEX(events[i]);
    lower[i] = boxes[b].min < cutOff;
    upper[i] = boxes[b].max > cutOff;
		      });

  return eventsPair(pack(events, lower),
		    pack(events, upper));
}

// n is the number of events (i.e. twice the number of triangles)
treeNode* generateNode(Boxes &boxes,
		       Events events,
		       BoundingBox B, 
		       size_t maxDepth) {
  timer t;
  index_t n = events[0].size();
  //cout << "event : " << events[0][0].v << ", " << IS_START(events[0][0]) << ", " << GET_INDEX(events[0][0]) << endl;
  //cout << "event : " << events[1][0].v << ", " << IS_START(events[1][0]) << ", " << GET_INDEX(events[1][0]) << endl;
  //cout << "n=" << n << " maxDepth=" << maxDepth << endl;
  if (n <= 2 || maxDepth == 0) 
    return treeNode::newNode(std::move(events), n, B);
  
  // loop over dimensions and find the best cut across all of them
  cutInfo cuts[3];
  parallel_for(0, 3, [&] (size_t d) {
    cuts[d] = bestCut(events[d], B[d], B[(d+1)%3], B[(d+2)%3]);
    //cout << "cost = " << cuts[d].cost << endl;
		     }, 1);

  //t.next("best cut");
  int cutDim = 0;
  for (int d = 1; d < 3; d++) 
    if (cuts[d].cost < cuts[cutDim].cost) cutDim = d;

  sequence<range>& cutDimRanges = boxes[cutDim];
  float cutOff = cuts[cutDim].cutOff;
  float area = boxSurfaceArea(B);
  float bestCost = CT + CL * cuts[cutDim].cost/area;
  float origCost = (float) (n/2);

  //cout << B << endl;
  //cout << bestCost << ", " << origCost << ", " << cutDim << ", " << area << ", " << CT << ", " << CL << ", " << cuts[cutDim].cost << endl;
  //abort();
  // quit recursion early if best cut is not very good
  if (bestCost >= origCost || 
      cuts[cutDim].numLeft + cuts[cutDim].numRight > maxExpand * n/2)
    return treeNode::newNode(std::move(events), n, B);

  // declare structures for recursive calls
  BoundingBox BBL;
  for (int i=0; i < 3; i++) BBL[i] = B[i];
  BBL[cutDim] = range(BBL[cutDim].min, cutOff);
  array<sequence<event>,3> leftEvents;
  index_t nl;

  BoundingBox BBR;
  for (int i=0; i < 3; i++) BBR[i] = B[i];
  BBR[cutDim] = range(cutOff, BBR[cutDim].max);
  array<sequence<event>,3> rightEvents;
  index_t nr;

  // now split each event array to the two sides
  eventsPair X[3];
  parallel_for (0, 3, [&] (size_t d) {
			//for (int d = 0; d < 3; d++) {
    X[d] = splitEvents(cutDimRanges, events[d], cutOff);
    },1);
  //t.next("split");
  
  for (int d = 0; d < 3; d++) {
    leftEvents[d] = std::move(X[d].first);
    rightEvents[d] = std::move(X[d].second);
    if (d == 0) {
      nl = X[d].first.size();
      nr = X[d].second.size();
    } else if (X[d].first.size() != nl || X[d].second.size() != nr) {
      cout << "kdTree: mismatched lengths, something wrong" << endl;
      abort();
    }
  }
  //t.next("left right");

  // free old events and make recursive calls
  for (int i=0; i < 3; i++) events[i] = sequence<event>();
  //t.next("clear");
  treeNode *L;
  treeNode *R;
  //abort();
  par_do([&] () {L = generateNode(boxes, std::move(leftEvents),
				  BBL, maxDepth-1);},
         [&] () {R = generateNode(boxes, std::move(rightEvents),
				  BBR, maxDepth-1);});

  return new treeNode(L, R, cutDim, cutOff, B);
  }

size_t tcount = 0;
size_t ccount = 0;

// Given an a ray, a bounding box, and a sequence of triangles, returns the 
// index of the first triangle the ray intersects inside the box.
// The triangles are given by n indices I into the triangle array Tri.
// -1 is returned if there is no intersection
index_t findRay(ray_t r,
	       sequence<index_t> const &I, 
	       triangles<point> const &Tri,
	       BoundingBox B) {
  index_t n = I.size();
  //cout << "Isize: " << I.size() << endl;
  if (STATS) { tcount += n; ccount += 1;}
  coord tMin = std::numeric_limits<double>::max();
  index_t k = -1;
  for (size_t i = 0; i < n; i++) {
    index_t j = I[i];
    point m[3] = {Tri.P[Tri.T[j][0]],  Tri.P[Tri.T[j][1]],  Tri.P[Tri.T[j][2]]};
    coord t = rayTriangleIntersect(r, m);
    if (t > 0.0 && t < tMin && inBox(r.o + r.d*t, B)) {
      tMin = t;
      k = j;
    }
  }
  return k;
}

// Given a ray and a tree node find the index of the first triangle the 
// ray intersects inside the box represented by that node.
// -1 is returned if there is no intersection
index_t findRay(ray_t r, treeNode* TN, triangs const &Tri) {
  //cout << "TN->n=" << TN->n << endl;
  if (TN->isLeaf()) 
    return findRay(r, TN->triangleIndices, Tri, TN->box);
  point o = r.o;
  vect d = r.d;

  coord oo[3] = {o.x,o.y,o.z};
  coord dd[3] = {d.x,d.y,d.z};

  // intersect ray with splitting plane
  int k0 = TN->cutDim;
  int k1 = (k0 == 2) ? 0 : k0+1;
  int k2 = (k0 == 0) ? 2 : k0-1;
  point2d<coord> o_p(oo[k1], oo[k2]);
  vector2d<coord> d_p(dd[k1], dd[k2]);
  // does not yet deal with dd[k0] == 0
  coord scale = (TN->cutOff - oo[k0])/dd[k0];
  point2d<coord> p_i = o_p + d_p * scale;

  range rx = TN->box[k1];
  range ry = TN->box[k2];
  coord d_0 = dd[k0];
  //cout << k0 << ", " << k1 << ", " << k2 << ", " << o_p << ", " << d_p << ", " << scale << ", " << rx.min << ", " << rx.max << ", " << ry.min << ", " << ry.max << ", " << d_0 << endl;

  // decide which of the two child boxes the ray intersects
  enum {LEFT, RIGHT, BOTH};
  int recurseTo = LEFT;
  if      (p_i.x < rx.min) { if (d_p.x*d_0 > 0) recurseTo = RIGHT;}
  else if (p_i.x > rx.max) { if (d_p.x*d_0 < 0) recurseTo = RIGHT;}
  else if (p_i.y < ry.min) { if (d_p.y*d_0 > 0) recurseTo = RIGHT;}
  else if (p_i.y > ry.max) { if (d_p.y*d_0 < 0) recurseTo = RIGHT;}
  else recurseTo = BOTH;
  //cout << "recurse to: " << recurseTo << ", " << TN->left->n << ", " << TN->right->n << endl;

  if (recurseTo == RIGHT) return findRay(r, TN->right, Tri);
  else if (recurseTo == LEFT) return findRay(r, TN->left, Tri);
  else if (d_0 > 0) {
    index_t t = findRay(r, TN->left, Tri);
    if (t >= 0) return t;
    else return findRay(r, TN->right, Tri);
  } else {
    index_t t = findRay(r, TN->right, Tri);
    if (t >= 0) return t;
    else return findRay(r, TN->left, Tri);
  }
}

sequence<index_t> rayCast(triangles<point> const &Tri,
		       sequence<ray<point>> const &rays) {
  timer t;
  index_t numRays = rays.size();

  // Extract triangles into a separate array for each dimension with
  // the lower and upper bound for each triangle in that dimension.
  Boxes boxes;
  index_t n = Tri.T.size();
  for (int d = 0; d < 3; d++)
    boxes[d] = sequence<range>::uninitialized(n);
  //sequence<point> P = Tri.P;
  parallel_for (0, n, [&] (size_t i) {
    point p0 = Tri.P[Tri.T[i][0]];
    point p1 = Tri.P[Tri.T[i][1]];
    point p2 = Tri.P[Tri.T[i][2]];
    boxes[0][i] = fixRange(min(p0.x,min(p1.x,p2.x)),max(p0.x,max(p1.x,p2.x)));
    boxes[1][i] = fixRange(min(p0.y,min(p1.y,p2.y)),max(p0.y,max(p1.y,p2.y)));
    boxes[2][i] = fixRange(min(p0.z,min(p1.z,p2.z)),max(p0.z,max(p1.z,p2.z)));
  });

  // Loop over the dimensions creating an array of events for each
  // dimension, sorting each one, and extracting the bounding box
  // from the first and last elements in the sorted events in each dim.
  Events events;
  BoundingBox boundingBox;
  for (int d = 0; d < 3; d++) {
    events[d] = tabulate(2*n, [&] (size_t i) -> event {
      return ((i % 2 == 0) ?
	      event(boxes[d][i/2].min, i/2, START) :
	      event(boxes[d][i/2].max, i/2, END));});
    sort_inplace(events[d], [] (event a, event b) {return a.v < b.v;}); 
    //std::sort(events[d].begin(), events[d].end(), [] (event a, event b) {return a.v < b.v;}); 
    boundingBox[d] = range(events[d][0].v, events[d][2*n-1].v);
  }
  // cout << boundingBox << endl;
  t.next("generate and sort events");

  // build the tree
  size_t recursionDepth = min<size_t>(maxRecursionDepth, parlay::log2_up(n)-1);
  treeNode* R = generateNode(boxes, std::move(events), boundingBox, recursionDepth);
  t.next("build tree");

  if (STATS)
    cout << "Triangles across all leaves = " << R->n 
	 << " Leaves = " << R->leaves << endl;

  tcount = 0;
  ccount = 0;

  // get the intersections
  auto results = tabulate(numRays, [&] (size_t i) -> index_t {
			       //cout << rays[i].o << " :: " << rays[i].d <<  endl;
			       return  findRay(rays[i], R, Tri);});
  t.next("intersect rays");

  treeNode::delete_tree(R);
  t.next("delete tree");
			  
  if (STATS)
    cout << "tcount=" << tcount << " ccount=" << ccount << endl;

  if (CHECK) {
    int nr = 10;
    auto indx = tabulate(n, [&] (size_t i) -> index_t {return i;});
    for (int i= 0; i < nr; i++) {
      //cout << results[i] << endl;
      if (findRay(rays[i], indx, Tri, boundingBox) != results[i]) {
	cout << "bad intersect in checking ray intersection" << endl;
	abort();
      }
    }
    t.next("check");
  }

  return results;
}
