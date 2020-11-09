//Timothy Chan 12/05
//approximate nearest neighbors: the SSS method (static version)

#include<iostream>
#include<cstdlib>
#include<math.h>
#include<limits>
#include<cfloat>
#include <algorithm>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/geometry.h"
#include "nearestNeighbors/octTree/oct_tree.h"
#define sq(x) (((float) (x))* ((float) (x)))
#define MAX (1<<29)

using o_tree = oct_tree<vtx>;

typedef int* Point; 
int d, shift; 
float eps, r, r_sq; 
Point ans, q1, q2; 
int key_bits = 64;
report_stats = true;

inline int less_msb(int x, int y) { return x < y && x < (x^y); }

//compares the interleaved bits of shifted points p, q without explicitly interleaving them
int cmp_shuffle(Point *p, Point *q){
	int j, k, x, y;
	for (j = k = x = 0; k<d; k++){
		if (less_msb(x, y = ((*p)[k]+shift)^(*q[k]+shift))){
			j=k; 
			x=y;
		}
	}
	return (*p)[j] - (*q)[k];
}

Point P[] convert(parlay::sequence(vtx*) &v, size_t n, int d){
	//prelims for rounding each point to an integer: 
	// 1) find the smallest point in each dimension
	// 2) find the largest gap between min and max over all dimensions
	using box = typename o_tree::box;
	box b = o_tree::get_box(V);
	double Delta = 0;
	for (int i = 0; i < d; i++) 
	  Delta = std::max(Delta, b.second[i] - b.first[i]); 
		min_point = b.first;
	//round each point to an integer with key_bits bit
	int bits = key_bits/d;
	uint maxval = (((size_t) 1) << bits) - 1; 
	uint ip[d];
    for (i=0; i<n; i++){
    	for (int j = 0; j < d; j++) 
      		ip[j] = floor(maxval * (p[j] - min_point[j])/delta); 
    }
    //now, convert each vertex into a point in Chan's implementation
    Point *P;
	P = new Point[n]; //idk if n can be of type size_t
	for(i = 0; i < n; i++){
		P[i] = new int[d];
		if (d == 2){
			P[i][0] = (i -> pt).x;
			P[i][1] = (i -> pt).y;
		} else if (d ==3){
			P[i][0] = (i -> pt).x;
			P[i][1] = (i -> pt).y;
			P[i][1] = (i -> pt).z;
		}
	}
	return P;
}

//sort the points based on their interleave ordering
void SSS_preprocess(Point P[], size_t n, int d){ //TODO check I'm not overloading d
	shift = (int) (drand48()*MAX);	
	q1 = new int[d];
	q2 = new int[d];
	qsort((void *) P, n, sizeof(Point),
		(int (*)(const void *, const void *)) cmp_shuffle);
}


void check_dist(Point p, Point q){
	int j; float z; 
	for (j=0, z=0; j<d; j++) z+= sq(p[j]-q[j]);
	if (z < r_sq) {
		r_sq = z; r = sqrt(z); ans = p;
		for (j=0; j<d; j++) {
			q1[j] = (q[j]<r) ? (q[j] - (int)ceil(r)): MAX;
		}
	}
}

//find the distance between query point q and the box that p is in
float dist_sq_to_box(Point q, Point p1, Point p2){
	int i, j, x, y; float z;
	for (j=x=0; j<d; j++)
		if (less_msb(x, y=(p1[j]+shift)^(p2[j]+shift))) x=y;
	frexp(x, &i);
	for (j=0, z=0; j<d; j++){
		x = ((p1[j]+shift)>>i)<<i; y = x+(1<<i);
		if (q[j]+shift < x) z += sq(q[j]+shift-x);
		else if (q[j] + shift > y) z += sq(q[j]+shift - y);
	}
	return z;
}

void SSS_query0(Point P[], size_t n, Point q){
	if (n==0) return;
	check_dist(P[n/2], q);
	if (n == 1 || dist_sq_to_box(q, P[0], P[n-1])*sq(1+eps) > r_sq) return;
	if (cmp_shuffle(&q, &P[n/2]) > 0) {
		SSS_query0(P, n/2, q);
		if (cmp_shuffle(&q2, &P[n/2]) < 0) SSS_query0(P + n/2+1, n-n/2-1, q);
	} else{
		SSS_query0(P + n/2+1, n-n/2+1, q);
		if (cmp_shuffle(&q1, &P[n/2]) < 0) SSS_query0(P, n/2, q);
	}
}

Point SSS_query(Point P[], int n, Point q){
	r_sq = FLT_MAX;
	SSS_query0(P, n, q);
	return ans;
}

template<class vtx>
void ANN(parlay::sequence(vtx*) &v, float eps){
	timer t("ANN", report_stats);{
		t.next("preprocess");
		size_t n = v.size(); 
		int dims = (v[0]->pt.dimension());
		srand48(12121+n+n+dims); //may not work with n as size_t
		Point P[] = convert(v, n, dims)
		SSS_preprocess(P, n, dims); 
		t.next("find all");
		parlay::parallel_for(0, n, [&] (size_t i)){
			SSS_query(P, n, i);
		}
	};
}

//questions about compiling using parlay:
// 1) what exactly does the makefile do here? do I just need what I have in it currently?
// 2) in testInputs, what do bnchmark and benchmark mean? why is there an identical copy in nearestNeighbors/bench?
// 3) in the octTree folder, what is the huge file just called neighbors?
// 4) it seems like some of the stuff in neighborsCheck and neighborsTime is needed here, what exactly?