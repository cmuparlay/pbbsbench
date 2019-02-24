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
#include "../pbbslib/stlalgs.h"
#include "../pbbslib/random_shuffle.h"

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
  
  pbbs::sequence<edge<intT>> E(M.nonZeros);
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
  return edgeArray<intT>(E,M.numRows,M.numCols);
}

edgeArray<intT> randomShuffle(edgeArray<intT> const &A) {
  auto E =  pbbs::random_shuffle(A.E);
  return edgeArray<intT>(E, A.numRows, A.numCols);
}

template <class intT>
edgeArray<intT> remDuplicates(edgeArray<intT> const &A) {
  auto lessE = [&] (edge<intT> a, edge<intT> b) {
    return (a.u < b.u) || ((a.u == b.u) && (a.v < b.v));};
  pbbs::sequence<edge<intT>> E =
    pbbs::remove_duplicates_ordered(A.E, lessE);
  return edgeArray<intT>(E, A.numRows, A.numCols);
}

template <class intT>
struct nEQF {bool operator() (edge<intT> e) const {return (e.u != e.v);}};

template <class intT>
edgeArray<intT> makeSymmetric(edgeArray<intT> const &A) {
  sequence<edge<intT>> EF = pbbs::filter(A.E, [&] (edge<intT> e) {
      return e.u != e.v;});
  auto FE = pbbs::delayed_seq<edge<intT>>(EF.size(), [&] (size_t i) {
      return edge<intT>(EF[i].v, EF[i].u);});
  return remDuplicates(edgeArray<intT>(pbbs::append(EF, FE),
				       A.numRows, A.numCols));
}

template <class intT>
struct getuF {intT operator() (edge<intT> e) const {return e.u;} };

template <class intT>
struct getuEdge {intT operator() (wghEdge<intT> e) const {return e.u;} };

template <class intT>
graph<intT> graphFromEdges(edgeArray<intT> const &EA, bool makeSym) {
  edgeArray<intT> A;
  if (makeSym) A = makeSymmetric<intT>(EA);
  else A = EA;
  intT m = A.nonZeros;
  intT n = max<intT>(A.numCols,A.numRows);
  sequence<size_t> offsets = intSort::iSort(A.E.begin(),m,n,getuF<intT>());
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
  pbbs::sequence<edge<intT>> E(nonZeros);
  intT k = 0;
  auto degrees = pbbs::delayed_seq<size_t>(numRows, [&] (size_t i) {
      return V[i].degree;});
  auto offsets = pbbs::scan(degrees, pbbs::addm<size_t>());
  parallel_for(0, numRows, [&] (size_t j) {
      size_t off = offsets.first[j];
      for (intT i = 0; i < V[j].degree; i++)
	E[off+i] = edge<intT>(j,V[j].Neighbors[i]);
    });
  return edgeArray<intT>(E,numRows,numRows);
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
graph<intT> graphReorder(graph<intT> Gr, sequence<intT> I = sequence<intT>(0)) {
  intT n = Gr.n;
  intT m = Gr.m;
  bool noI = (I.size()==0);
  if (noI) I = pbbs::random_permutation<intT>(n);
  vertex<intT> *V = newA(vertex<intT>,Gr.n);
  parallel_for(0, n, [&] (size_t i) {
      V[I[i]] = Gr.V[i];});
  parallel_for (0, Gr.n, [&] (size_t i) {
    for (intT j=0; j < V[i].degree; j++) {
      V[i].Neighbors[j] = I[V[i].Neighbors[j]];
    }
    std::sort(V[i].Neighbors,V[i].Neighbors+V[i].degree);
    });
  free(Gr.V);
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

template <typename intT>
void writeGraphFromAdj(graph<intT> G, char* fname, bool adjArray, bool ordered) {
  if (!ordered) G = graphReorder<intT>(G);
  if (adjArray) writeGraphToFile<intT>(G, fname);
  else {
    auto B = edgesFromGraph<intT>(G);
    if (!ordered) B = randomShuffle(B);
    writeEdgeArrayToFile<intT>(B, fname);
  }
}

template <typename intT>
void writeGraphFromEdges(edgeArray<intT> const & EA, char* fname, bool adjArray, bool ordered) {
  graph<intT> G = graphFromEdges<intT>(EA, adjArray);
  writeGraphFromAdj(G, fname, adjArray, ordered);
}

// template <class intT>
// edgeArray<intT> edgesFromMtxFile(const char* fname) {
//   ifstream file (fname, ios::in);
//   char* line = newA(char,1000);
//   intT i,j = 0;
//   while (file.peek() == '%') {
//     j++;
//     file.getline(line,1000);
//   }
//   intT numRows, numCols, nonZeros;
//   file >> numRows >> numCols >> nonZeros;
//   //cout << j << "," << numRows << "," << numCols << "," << nonZeros << endl;
//   edge<intT> *E = newA(edge<intT>,nonZeros);
//   double toss;
//   for (i=0, j=0; i < nonZeros; i++) {
//     file >> E[j].u >> E[j].v >> toss;
//     E[j].u--;
//     E[j].v--;
//     if (toss != 0.0) j++;
//   }
//   nonZeros = j;
//   //cout << "nonzeros = " << nonZeros << endl;
//   file.close();  
//   return edgeArray<intT>(E,numRows,numCols,nonZeros);
// }

