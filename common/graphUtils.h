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
#include "../pbbslib/parallel.h"
#include "../pbbslib/quicksort.h"
#include "../pbbslib/stlalgs.h"
#include "../pbbslib/random_shuffle.h"
#include "../pbbslib/integer_sort.h"

using namespace std;

namespace dataGen {

#define HASH_MAX_INT ((unsigned) 1 << 31)

  //#define HASH_MAX_LONG ((unsigned long) 1 << 63)

  template <class T> T hash(intT i);
  
  template <>
  inline intT hash<intT>(intT i) {
    return pbbs::hash32(i) & (HASH_MAX_INT-1);}

  template <>
  inline uintT hash<uintT>(intT i) {
    return pbbs::hash32(i);}

  template <>
  inline double hash<double>(intT i) {
    return ((double) hash<intT>(i)/((double) HASH_MAX_INT));}

};

template <class intT>
wghEdgeArray<intT> addRandWeights(edgeArray<intT> const &G) {
  pbbs::random r(257621);
  intT m = G.nonZeros;
  intT n = G.numRows;
  sequence<wghEdge<intT>> E(m, [&] (size_t i) {
      return wghEdge<intT>(G.E[i].u, G.E[i].v, dataGen::hash<double>(i));});
  return wghEdgeArray<intT>(std::move(E), n);
}

// template <class intT>
// edgeArray<intT> edgesFromSparse(sparseRowMajor<double,intT> M) {
//   pbbs::sequence<edge<intT>> E(M.nonZeros);
//   intT k = 0;
//   for (intT i=0; i < M.numRows; i++) {
//     for (intT j=M.Starts[i]; j < M.Starts[i+1]; j++) {
//       if (M.Values[j] != 0.0) {
// 	E[k].u = i;
// 	E[k].v = M.ColIds[j];
// 	k++;
//       }
//     }
//   }
//   intT nonZeros = k;
//   return edgeArray<intT>(std::move(E), M.numRows, M.numCols);
// }

edgeArray<intT> randomShuffle(edgeArray<intT> const &A) {
  auto E =  pbbs::random_shuffle(A.E);
  return edgeArray<intT>(std::move(E), A.numRows, A.numCols);
}

template <class intT>
edgeArray<intT> remDuplicates(edgeArray<intT> const &A) {
  auto lessE = [&] (edge<intT> a, edge<intT> b) {
    return (a.u < b.u) || ((a.u == b.u) && (a.v < b.v));};
  pbbs::sequence<edge<intT>> E =
    pbbs::remove_duplicates_ordered(A.E, lessE);
  return edgeArray<intT>(std::move(E), A.numRows, A.numCols);
}

template <class intT>
edgeArray<intT> makeSymmetric(edgeArray<intT> const &A) {
  pbbs::sequence<edge<intT>> EF = pbbs::filter(A.E, [&] (edge<intT> e) {
      return e.u != e.v;});
  auto FE = pbbs::delayed_seq<edge<intT>>(EF.size(), [&] (size_t i) {
      return edge<intT>(EF[i].v, EF[i].u);});
  return remDuplicates(edgeArray<intT>(pbbs::append(EF, FE),
				       A.numRows, A.numCols));
}

template <class intT>
struct getuEdge {intT operator() (wghEdge<intT> e) const {return e.u;} };

template <class intT>
graph<intT> graphFromEdges(edgeArray<intT> const &EA, bool makeSym) {
  edgeArray<intT> SA;
  if (makeSym) SA = makeSymmetric<intT>(EA);
  edgeArray<intT> const &A = (makeSym) ? SA : EA;

  intT m = A.nonZeros;
  intT n = std::max(A.numCols, A.numRows);

  pbbs::sequence<size_t> counts;
  pbbs::sequence<intT> offsets;
  pbbs::sequence<edge<intT>> E;
  size_t nn;
  auto getu = [&] (edge<intT> e) {return e.u;};
  std::tie(E, counts) = pbbs::integer_sort_with_counts(A.E, getu, n);
  std::tie(offsets,nn) = pbbs::scan(pbbs::delayed_seq<intT>(n+1, [&] (size_t i) {
	return (i == n) ? 0 : counts[i];}), pbbs::addm<intT>());

  return graph<intT>(std::move(offsets),
		     sequence<intT>(m, [&] (size_t i) {return E[i].v;}),
		     n);
}

template <class intT>
wghGraph<intT> wghGraphFromEdges(wghEdgeArray<intT> const &A) {
  intT n = A.n;
  intT m = A.m;

  pbbs::sequence<size_t> counts;
  pbbs::sequence<intT> offsets;
  pbbs::sequence<wghEdge<intT>> E;
  size_t nn;
  auto getu = [&] (wghEdge<intT> e) {return e.u;};
  std::tie(E, counts) = pbbs::integer_sort_with_counts(A.E, getu, n);
  std::tie(offsets,nn) = pbbs::scan(pbbs::delayed_seq<intT>(n+1, [&] (size_t i) {
	return (i == n) ? 0 : counts[i];}), pbbs::addm<intT>());

  return wghGraph<intT>(std::move(offsets),
			sequence<intT>(m, [&] (size_t i) {return E[i].v;}),
			sequence<intT>(m, [&] (size_t i) {
			    return (intT) 1000000000 * E[i].weight;}),
			n);
}

template <class intT>
edgeArray<intT> edgesFromGraph(graph<intT> const &G) {
  intT numRows = G.numVertices();
  intT nonZeros = G.numEdges();

  // flatten
  pbbs::sequence<edge<intT>> E(nonZeros);
  parallel_for(0, numRows, [&] (size_t j) {
      size_t off = G.get_offsets()[j];
      vertex<intT> v = G[j];
      for (intT i = 0; i < v.degree; i++)
	E[off+i] = edge<intT>(j, v.Neighbors[i]);
    });
  return edgeArray<intT>(std::move(E), numRows, numRows);
}

// template <class eType, class intT>
// sparseRowMajor<eType,intT> sparseFromGraph(graph<intT> G) {
//   intT numRows = G.n;
//   intT nonZeros = G.m;
//   vertex<intT>* V = G.V;
//   intT *Starts = newA(intT,numRows+1);
//   intT *ColIds = newA(intT,nonZeros);
//   intT start = 0;
//   for (intT i = 0; i < numRows; i++) {
//     Starts[i] = start;
//     start += V[i].degree;
//   }
//   Starts[numRows] = start;
//   parallel_for (0, numRows, [&] (size_t j) {
//     for (intT i = 0; i < (Starts[j+1] - Starts[j]); i++) {
//       ColIds[Starts[j]+i] = V[j].Neighbors[i];
//     }});
//   return sparseRowMajor<eType,intT>(numRows,numRows,nonZeros,Starts,ColIds,NULL);
// }

// offset for start of each vertex if flattening the edge listd
template <class intT>
sequence<intT> getOffsets(sequence<vertex<intT>> const &V) {
  size_t n = V.size();
  auto degrees = pbbs::delayed_seq<intT>(n+1, [&] (size_t i) {
      return (i == n) ? 0 : V[i].degree;});
  return pbbs::scan(degrees, pbbs::addm<intT>()).first;
}

// if I is NULL then it randomly reorders
template <class intT>
graph<intT> graphReorder(graph<intT> const &Gr,
			 pbbs::sequence<intT> const &I = pbbs::sequence<intT>(0)) {
  intT n = Gr.numVertices();
  intT m = Gr.numEdges();

  bool noI = (I.size()==0);
  pbbs::sequence<intT> const &II = noI ? pbbs::random_permutation<intT>(n) : I;

  // now write vertices to new locations
  // inverse permutation
  pbbs::sequence<vertex<intT>> V(n);
  parallel_for (0, n, [&] (size_t i) {
      V[II[i]] = Gr[i];});
  pbbs::sequence<intT> offsets = getOffsets(V);
  pbbs::sequence<intT> E(m);
  parallel_for (0, n, [&] (size_t i) {
      size_t o = offsets[i];
      for (intT j=0; j < V[i].degree; j++) 
	E[o + j] = II[V[i].Neighbors[j]];
      std::sort(E.begin() + o, E.begin() + o + V[i].degree);
    }, 1000);
  return graph<intT>(std::move(offsets), std::move(E), n);
}

template <class intT>
int graphCheckConsistency(graph<intT> const &Gr) {
  size_t n = Gr.numVertices();
  size_t m = Gr.numEdges();
  size_t edgecount = pbbs::reduce(pbbs::delayed_seq<size_t>(n, [&] (size_t i) {
	return Gr[i].degree;}), pbbs::addm<size_t>());
  if (m != edgecount) {
    cout << "bad edge count in graphCheckConsistency: m = " 
	 << m << " sum of degrees = " << edgecount << endl;
    return 1;
  }
  size_t error_loc = pbbs::reduce(pbbs::delayed_seq<size_t>(n, [&] (size_t i) {
	for (intT j=0; j < Gr[i].degree; j++) 
	  if (Gr[i].Neighbors[j] >= n) return i;
	return n;
      }), pbbs::minm<size_t>());
  if (error_loc < n) {
    cout << "edge out of range in graphCheckConsistency: at i = " 
	 << error_loc << endl;
    return 1;
  }
}

// template <class intT>
// sparseRowMajor<double,intT> sparseFromCsrFile(const char* fname) {
//   FILE *f = fopen(fname,"r");
//   if (f == NULL) {
//     cout << "Trying to open nonexistant file: " << fname << endl;
//     abort();
//   }

//   intT numRows;  intT numCols;  intT nonZeros;
//   intT nc = fread(&numRows, sizeof(intT), 1, f);
//   nc = fread(&numCols, sizeof(intT), 1, f);
//   nc = fread(&nonZeros, sizeof(intT), 1, f); 

//   double *Values = (double *) malloc(sizeof(double)*nonZeros);
//   intT *ColIds = (intT *) malloc(sizeof(intT)*nonZeros);
//   intT *Starts = (intT *) malloc(sizeof(intT)*(1 + numRows));
//   Starts[numRows] = nonZeros;

//   size_t r;
//   r = fread(Values, sizeof(double), nonZeros, f);
//   r = fread(ColIds, sizeof(intT), nonZeros, f);
//   r = fread(Starts, sizeof(intT), numRows, f); 
//   fclose(f);
//   return sparseRowMajor<double,intT>(numRows,numCols,nonZeros,Starts,ColIds,Values);
// }

// The following two are used by the graph generators to write out in either format
// and either with reordering or not
template <typename intT>
void writeGraphFromAdj(graph<intT> const &G, char* fname, bool adjArray, bool ordered) {
  auto GR = (!ordered) ? graphReorder<intT>(G) : G;
  if (adjArray) writeGraphToFile<intT>(G, fname);
  else {
    auto B = edgesFromGraph<intT>(G);
    if (!ordered) B = randomShuffle(B);
    writeEdgeArrayToFile<intT>(B, fname);
  }
}

template <typename intT>
void writeGraphFromEdges(edgeArray<intT> const & EA, char* fname, bool adjArray, bool ordered) {
  graph<intT> const &G = graphFromEdges<intT>(EA, adjArray);
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

