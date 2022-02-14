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

#include <algorithm>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/geometry.h"


template<class fvec_point>
float distance(fvec_point* p, fvec_point* q, unsigned d){
	float dist = 0;
	for(int i=0; i<d; i++) dist += ((p->coordinates)[i] - (q->coordinates[i]))*((p->coordinates)[i] - (q->coordinates[i]));
	return dist;
}

// template<class fvec_point>
// float distance(fvec_point* p, fvec_point* q){
// 	float dist = 0;

// 	float diff0, diff1, diff2, diff3;

// 	int d = ((p->coordinates).size())/4;
// 	int i=0;
// 	while(i<d-3){
// 		diff0 = (p->coordinates)[i] - (q->coordinates)[i];
// 		diff1 = (p->coordinates)[i+1] - (q->coordinates)[i+1];
// 		diff2 = (p->coordinates)[i+2] - (q->coordinates)[i+2];
// 		diff3 = (p->coordinates)[i+3] - (q->coordinates)[i+3];
// 		dist += diff0 * diff0 + diff1 * diff1 + diff2 * diff2 + diff3 * diff3;
// 		i+=4;
// 	}
// 	float rem;
// 	for(int i; i<d; i++){
// 		rem = (p->coordinates)[i] - (q->coordinates)[i];
// 		dist += rem*rem;
// 	}
// 	return dist; 
// }