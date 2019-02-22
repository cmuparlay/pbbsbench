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

#pragma once

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include "graph.h"
#include "glue.h"
#include "../pbbslib/parallel.h"
#include "../pbbslib/quicksort.h"

using namespace std;

template <class intT>
wghEdgeArray<intT> addRandWeights(edgeArray<intT> G) {
  intT m = G.nonZeros;
  intT n = G.numRows;
  wghEdge<intT> *E = newA(wghEdge<intT>, m);
  for (intT i=0; i < m; i++) {
    E[i].u = G.E[i].u;
    E[i].v = G.E[i].v;
    E[i].weight = utils::hashInt(i);
  }
  return wghEdgeArray<intT>(E, n, m);
}

template <class intT>
edgeArray<intT> edgesFromSparse(sparseRowMajor<double,intT> M) {
  edge<intT> *E = newA(edge<intT>,M.nonZeros);
  intT k = 0;
  for (intT i=0; i < M.numRows; i++) {
    for (intT j=M.Starts[i]; j < M.Starts[i+1]; j++) {
      if (M.Values[j] != 0.0) {
	E[k].u = i;
	E[k].v = M.ColIds[j];
	k++;
      }
    }
  }
  intT nonZeros = k;
  return edgeArray<intT>(E,M.numRows,M.numCols,nonZeros);
}

template <class intT>
int cmpInt(intT v, intT b) {
  return (v > b) ? 1 : ((v == b) ? 0 : -1);}

template <class intT>
struct hashEdge {
  typedef edge<intT>* eType;
  typedef edge<intT>* kType;
  eType empty() {return NULL;}
  kType getKey(eType v) {return v;}
  intT hash(kType e) {
    return utils::hashInt(e->u) + utils::hashInt(100*e->v); }
  int cmp(kType a, kType b) {
    int c = cmpInt(a->u, b->u);
    return (c == 0) ? cmpInt(a->v,b->v) : c;
  }
  bool replaceQ(eType v, eType b) {return 0;}
};

template <class intT>
_seq<edge<intT>* > removeDuplicates(_seq<edge<intT> *>  S) {
  return removeDuplicates(S,hashEdge<intT>());}



template <class intT>
edgeArray<intT> remDuplicates(edgeArray<intT> A) {
  intT m = A.nonZeros;
  edge<intT> **EP = newA(edge<intT>*,m);
  parallel_for (0, m, [&] (size_t i) {EP[i] = A.E+i;});
   _seq<edge<intT> *> F = removeDuplicates(_seq<edge<intT> *>(EP,m));
   //_seq<edge<intT>* > F = removeDuplicates(_seq<edge<intT> *>(EP,m),hashEdge<intT>());
  free(EP);
  intT l = F.n;
  edge<intT> *E = newA(edge<intT>,m);
  parallel_for (0, l, [&] (size_t j) {E[j] = *F.A[j];});
  F.del();
  return edgeArray<intT>(E,A.numRows,A.numCols,l);
}

struct edgeCmp {
  bool operator() (edge<intT> e1, edge<intT> e2) const {
    return ((e1.u < e2.u) ? 1 : ((e1.u > e2.u) ? 0 : (e1.v < e2.v)));
  }
};


template <class intT>
edgeArray<intT> remDuplicates2(edgeArray<intT> A) {
  intT m = A.nonZeros;
  edge<intT> * E = newA(edge<intT>,m);
  parallel_for(0, m, [&] (size_t i) {
      E[i].u = A.E[i].u; E[i].v = A.E[i].v;});
  pbbs::quicksort(E,m,edgeCmp());
  //for(int i=0;i<m;i++)cout<<"("<<E[i].u<<","<<E[i].v<<") ";cout<<endl;

  intT* flags = newA(intT,m);
  flags[0] = 1;
  parallel_for(1, m, [&] (size_t i) {
    if((E[i].u != E[i-1].u) || (E[i].v != E[i-1].v)) flags[i] = 1;
    else flags[i] = 0;
    });

  intT mm = osequence::plusScan(flags,flags,m);

  //cout<<mm<<endl;
  edge<intT>* F = newA(edge<intT>,mm);
  F[mm-1] = E[m-1];
  parallel_for(0, m-1, [&] (size_t i) {
    if(flags[i] != flags[i+1]) F[flags[i]] = E[i];
    });
  free(flags);
  //for(int i=0;i<mm;i++)cout<<"("<<F[i].u<<","<<F[i].v<<") ";cout<<endl;

  return edgeArray<intT>(F,A.numRows,A.numCols,mm);
}


template <class intT>
struct nEQF {bool operator() (edge<intT> e) const {return (e.u != e.v);}};

template <class intT>
edgeArray<intT> makeSymmetric(edgeArray<intT> A) {
  intT m = A.nonZeros;
  edge<intT> *E = A.E;
  edge<intT> *F = newA(edge<intT>,2*m);
  intT mm = osequence::filter(E,F,m,nEQF<intT>());

  parallel_for (0, mm, [&] (size_t i) {
    F[i+mm].u = F[i].v;
    F[i+mm].v = F[i].u;
    });

  edgeArray<intT> R = remDuplicates2(edgeArray<intT>(F,A.numRows,A.numCols,2*mm));
  free(F);
  
  return R;
  //return edgeArray<intT>(F,A.numRows,A.numCols,2*mm);
}

template <class intT>
struct getuF {intT operator() (edge<intT> e) const {return e.u;} };

template <class intT>
struct getuEdge {intT operator() (wghEdge<intT> e) const {return e.u;} };

template <class intT>
graph<intT> graphFromEdges(edgeArray<intT> EA, bool makeSym) {
  edgeArray<intT> A;
  if (makeSym) A = makeSymmetric<intT>(EA);
  else {  // should have copy constructor
    edge<intT> *E = newA(edge<intT>,EA.nonZeros);
    parallel_for (0, EA.nonZeros, [&] (size_t i) {E[i] = EA.E[i];});
    A = edgeArray<intT>(E,EA.numRows,EA.numCols,EA.nonZeros);
  }
  intT m = A.nonZeros;
  intT n = max<intT>(A.numCols,A.numRows);
  sequence<size_t> offsets = intSort::iSort(A.E,m,n,getuF<intT>());
  intT *X = newA(intT,m);
  vertex<intT> *v = newA(vertex<intT>,n);
  parallel_for (0, n, [&] (size_t i) {
    intT o = offsets[i];
    intT l = ((i == n-1) ? m : offsets[i+1])-offsets[i];
    v[i].degree = l;
    v[i].Neighbors = X+o;
    for (intT j=0; j < l; j++) {
      v[i].Neighbors[j] = A.E[o+j].v;
    }
    });
  A.del();
  return graph<intT>(v,n,m,X);
}

template <class intT>
wghGraph<intT> wghGraphFromEdges(wghEdgeArray<intT> EA) {
  intT m = EA.m;
  intT n = EA.n;
  wghEdge<intT> *E = newA(wghEdge<intT>,m);
  parallel_for (0, m, [&] (size_t i) {E[i] = EA.E[i];});
  wghEdgeArray<intT> A = wghEdgeArray<intT>(E,n,m);
  sequence<size_t> offsets = intSort::iSort(A.E,m,n,getuEdge<intT>());
  intT *X = newA(intT,m);
  intT *weights = newA(intT,m);
  wghVertex<intT> *v = newA(wghVertex<intT>,n);
  parallel_for (0, n, [&] (size_t i) {
    intT o = offsets[i];
    intT l = ((i == n-1) ? m : offsets[i+1])-offsets[i];
    v[i].degree = l;
    v[i].Neighbors = X+o;
    v[i].nghWeights = weights+o;
    for (intT j=0; j < l; j++) {
      v[i].Neighbors[j] = A.E[o+j].v;
      v[i].nghWeights[j] = A.E[o+j].weight;
    }
    });
  A.del();
  return wghGraph<intT>(v,n,m,X,weights);
}

template <class intT>
edgeArray<intT> edgesFromGraph(graph<intT> G) {
  intT numRows = G.n;
  intT nonZeros = G.m;
  vertex<intT>* V = G.V;
  edge<intT> *E = newA(edge<intT>, nonZeros);
  intT k = 0;
  for (intT j=0; j < numRows; j++)
    for (intT i = 0; i < V[j].degree; i++)
      E[k++] = edge<intT>(j,V[j].Neighbors[i]);
  return edgeArray<intT>(E,numRows,numRows,nonZeros);
}

template <class eType, class intT>
sparseRowMajor<eType,intT> sparseFromGraph(graph<intT> G) {
  intT numRows = G.n;
  intT nonZeros = G.m;
  vertex<intT>* V = G.V;
  intT *Starts = newA(intT,numRows+1);
  intT *ColIds = newA(intT,nonZeros);
  intT start = 0;
  for (intT i = 0; i < numRows; i++) {
    Starts[i] = start;
    start += V[i].degree;
  }
  Starts[numRows] = start;
  parallel_for (0, numRows, [&] (size_t j) {
    for (intT i = 0; i < (Starts[j+1] - Starts[j]); i++) {
      ColIds[Starts[j]+i] = V[j].Neighbors[i];
    }});
  return sparseRowMajor<eType,intT>(numRows,numRows,nonZeros,Starts,ColIds,NULL);
}

// if I is NULL then it randomly reorders
template <class intT>
graph<intT> graphReorder(graph<intT> Gr, intT* I) {
  intT n = Gr.n;
  intT m = Gr.m;
  bool noI = (I==NULL);
  if (noI) {
    I = newA(intT,Gr.n);
    parallel_for (0, Gr.n, [&] (size_t i) {I[i] = i;});
    std::random_shuffle(I,I+Gr.n);
  }
  vertex<intT> *V = newA(vertex<intT>,Gr.n);
  for (intT i=0; i < Gr.n; i++) V[I[i]] = Gr.V[i];
  for (intT i=0; i < Gr.n; i++) {
    for (intT j=0; j < V[i].degree; j++) {
      V[i].Neighbors[j] = I[V[i].Neighbors[j]];
    }
    std::sort(V[i].Neighbors,V[i].Neighbors+V[i].degree);
  }
  free(Gr.V);
  if (noI) free(I);
  return graph<intT>(V,n,m,Gr.allocatedInplace);
}

template <class intT>
int graphCheckConsistency(graph<intT> Gr) {
  vertex<intT> *V = Gr.V;
  intT edgecount = 0;
  for (intT i=0; i < Gr.n; i++) {
    edgecount += V[i].degree;
    for (intT j=0; j < V[i].degree; j++) {
      intT ngh = V[i].Neighbors[j];
      utils::myAssert(ngh >= 0 && ngh < Gr.n,
		      "graphCheckConsistency: bad edge");
    }
  }
  if (Gr.m != edgecount) {
    cout << "bad edge count in graphCheckConsistency: m = " 
	 << Gr.m << " sum of degrees = " << edgecount << endl;
    abort();
  }
  return 0;
}

template <class intT>
sparseRowMajor<double,intT> sparseFromCsrFile(const char* fname) {
  FILE *f = fopen(fname,"r");
  if (f == NULL) {
    cout << "Trying to open nonexistant file: " << fname << endl;
    abort();
  }

  intT numRows;  intT numCols;  intT nonZeros;
  intT nc = fread(&numRows, sizeof(intT), 1, f);
  nc = fread(&numCols, sizeof(intT), 1, f);
  nc = fread(&nonZeros, sizeof(intT), 1, f); 

  double *Values = (double *) malloc(sizeof(double)*nonZeros);
  intT *ColIds = (intT *) malloc(sizeof(intT)*nonZeros);
  intT *Starts = (intT *) malloc(sizeof(intT)*(1 + numRows));
  Starts[numRows] = nonZeros;

  size_t r;
  r = fread(Values, sizeof(double), nonZeros, f);
  r = fread(ColIds, sizeof(intT), nonZeros, f);
  r = fread(Starts, sizeof(intT), numRows, f); 
  fclose(f);
  return sparseRowMajor<double,intT>(numRows,numCols,nonZeros,Starts,ColIds,Values);
}

template <class intT>
edgeArray<intT> edgesFromMtxFile(const char* fname) {
  ifstream file (fname, ios::in);
  char* line = newA(char,1000);
  intT i,j = 0;
  while (file.peek() == '%') {
    j++;
    file.getline(line,1000);
  }
  intT numRows, numCols, nonZeros;
  file >> numRows >> numCols >> nonZeros;
  //cout << j << "," << numRows << "," << numCols << "," << nonZeros << endl;
  edge<intT> *E = newA(edge<intT>,nonZeros);
  double toss;
  for (i=0, j=0; i < nonZeros; i++) {
    file >> E[j].u >> E[j].v >> toss;
    E[j].u--;
    E[j].v--;
    if (toss != 0.0) j++;
  }
  nonZeros = j;
  //cout << "nonzeros = " << nonZeros << endl;
  file.close();  
  return edgeArray<intT>(E,numRows,numCols,nonZeros);
}

