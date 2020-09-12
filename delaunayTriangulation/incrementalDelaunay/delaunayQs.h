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

// Code used in incremental Delaunay triangulation and incremental Delaunay refinement.
// For each vertex being inserted in parallel we store a stack of all the vertices and 
// simplices in its "tent" (i.e. the region that will be wiped out by the insersion).

#include "topology.h"

// A simplified version of STL vectors.   Fully upward compatible, but seems to be
// significantly faster when used in parallel.
template <class ET>
struct myVector {
  ET* vals;
  intT max, top;
  myVector() {
    top = 0;
    max = 50;
    vals = newA(ET, max);
  }
  void push_back(ET v) {
    if (top == max) {
      ET *x = newA(ET, 2*max);
      for (int j=0; j<top; j++) x[j] = vals[j];
      free(vals);
      vals = x;
      max = 2*max;
    }
    vals[top++] = v;
  }
  void clear() {top = 0;}
  intT size() { return top;}
  ET operator [] (intT i) {return vals[i];}
  void del() {free(vals);}
  ~myVector() {free(vals);}
};

// Holds vertex and simplex queues used to store the cavity created 
// while searching from a vertex between when it is initially searched 
// and later checked to see if all corners are reserved.
struct Qs {
  myVector<vertex*> vertexQ;
  myVector<simplex> simplexQ;
};

// A Queue for each vertex to be inserted in parallel to keep info
// between the reserve and commit phases.
// The size n puts a limit on how many can be run in parallel
struct vertexQs {
  Qs **qs;
  int n;
  vertexQs(int _n) {
    n = _n;
    qs = newA(Qs*, n);
    for (intT i=0; i < n; i++) qs[i] = new Qs;
  }
  void del() { 
    for (int i=0; i < n; i++) delete qs[i];
    free(qs); 
  }
  Qs* operator [] (intT i) {return qs[i];}
};
