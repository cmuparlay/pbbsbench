// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2010 Guy Blelloch and the PBBS team
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
#include <stdint.h>
#include <cstring>
#include "../pbbslib/parallel.h"
#include "IO.h"

#include <sys/mman.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

using namespace benchIO;

template <class intT>
int xToStringLen(edge<intT> a) {
  return xToStringLen(a.u) + xToStringLen(a.v) + 1;
}

template <class intT>
void xToString(char* s, edge<intT> a) {
  int l = xToStringLen(a.u);
  xToString(s, a.u);
  s[l] = ' ';
  xToString(s+l+1, a.v);
}

template <class intT>
int xToStringLen(wghEdge<intT> a) {
  return xToStringLen(a.u) + xToStringLen(a.v) + xToStringLen(a.weight) + 2;
}

template <class intT>
void xToString(char* s, wghEdge<intT> a) {
  int lu = xToStringLen(a.u);
  int lv = xToStringLen(a.v);
  xToString(s, a.u);
  s[lu] = ' ';
  xToString(s+lu+1, a.v);
  s[lu+lv+1] = ' ';
  xToString(s+lu+lv+2, a.weight);
}

namespace benchIO {
  using namespace std;

  string AdjGraphHeader = "AdjacencyGraph";
  string EdgeArrayHeader = "EdgeArray";
  string WghEdgeArrayHeader = "WeightedEdgeArray";
  string WghAdjGraphHeader = "WeightedAdjacencyGraph";

  template <class intT>
  int writeGraphToFile(graph<intT> const &G, char* fname) {
    intT m = G.numEdges();
    intT n = G.numVertices();
    intT totalLen = 2 + n + m;
    pbbs::sequence<intT> Out(totalLen);
    Out[0] = n;
    Out[1] = m;

    // write offsets to Out[2,..,2+n)
    auto offsets = G.get_offsets();
    parallel_for (0, n, [&] (size_t i) {
	Out[i+2] = offsets[i];});

    // write out edges to Out[2+n,..,2+n+m)
    parallel_for(0, n, [&] (size_t i) {
	size_t o = offsets[i] + 2 + n;
	vertex<intT> v = G[i];
	for (intT j = 0; j < v.degree; j++)
	  Out[o + j] = v.Neighbors[j];
      });
    int r = writeArrayToFile(AdjGraphHeader, Out.begin(), totalLen, fname);
    return r;
  }

  template <class intT>
  int writeWghGraphToFile(wghGraph<intT> G, char* fname) {
    intT m = G.m;
    intT n = G.n;
    intT totalLen = 2 + n + m*2;
    pbbs::sequence<intT> Out(totalLen);
    Out[0] = n;
    Out[1] = m;

    // write offsets to Out[2,..,2+n)
    auto offsets = G.get_offsets();
    parallel_for (0, n, [&] (size_t i) {
	Out[i+2] = offsets[i];});

    // write out edges to Out[2+n,..,2+n+m)
    // and weights to Out[2+n+m,..,2+n+2*m)
    parallel_for(0, n, [&] (size_t i) {
	size_t o = offsets[i] + 2 + n;
	wghVertex<intT> v = G[i];
	for (intT j = 0; j < v.degree; j++) {
	  Out[o + j] = v.Neighbors[j];
	  Out[o + m + j] = v.nghWeights[j]; }
      });
    int r = writeArrayToFile(WghAdjGraphHeader, Out.begin(), totalLen, fname);
    return r;
  }

  template <class intT>
  int writeEdgeArrayToFile(edgeArray<intT> const &EA, char* fname) {
    intT m = EA.nonZeros;
    int r = writeArrayToFile(EdgeArrayHeader, EA.E.begin(), m, fname);
    return r;
  }

  template <class intT>
  int writeWghEdgeArrayToFile(wghEdgeArray<intT> const &EA, char* fname) {
    uintT m = EA.m;
    int r = writeArrayToFile(WghEdgeArrayHeader, EA.E.begin(), m, fname);
    return r;
  }

  template <class intT>
  edgeArray<intT> readEdgeArrayFromFile(char* fname) {
    pbbs::sequence<char> S = readStringFromFile(fname);
    pbbs::sequence<char*> W = stringToWords(S);
    if (W[0] != EdgeArrayHeader) {
      cout << "Bad input file" << endl;
      abort();
    }
    long n = (W.size()-1)/2;
    pbbs::sequence<edge<intT>> E(n, [&] (long i) {
	return edge<intT>(atol(W[2*i + 1]),
			  atol(W[2*i + 2]));});

    auto mon = pbbs::make_monoid([&] (edge<intT> a, edge<intT> b) {
	return edge<intT>(std::max(a.u, b.u), std::max(a.v, b.v));},
      edge<intT>(0,0));
    auto r = pbbs::reduce(E, mon);

    intT maxrc = std::max(r.u, r.v) + 1;
    return edgeArray<intT>(E, maxrc, maxrc);
  }

  template <class intT>
  wghEdgeArray<intT> readWghEdgeArrayFromFile(char* fname) {
    pbbs::sequence<char> S = readStringFromFile(fname);
    pbbs::sequence<char*> W = stringToWords(S);
    if (W[0] != WghEdgeArrayHeader) {
      cout << "Bad input file" << endl;
      abort();
    }
    long n = (W.size()-1)/3;
    pbbs::sequence<wghEdge<intT>> E(n, [&] (size_t i) {
	return wghEdge<intT>(atol(W[3*i + 1]),
			     atol(W[3*i + 2]),
			     atof(W[3*i + 3]));});

    auto mon = pbbs::make_monoid([&] (wghEdge<intT> a, wghEdge<intT> b) {
	return wghEdge<int>(std::max(a.u, b.u), std::max(a.v, b.v), 0);},
      wghEdge<intT>(0,0,0));
    auto r = pbbs::reduce(E, mon);

    return wghEdgeArray<intT>(E, max<intT>(r.u, r.v) + 1);
  }

  template <class intT>
  graph<intT> readGraphFromFile(char* fname) {
    pbbs::sequence<char> S = readStringFromFile(fname);
    pbbs::sequence<char*> W = stringToWords(S);
    if (W[0] != AdjGraphHeader) {
      cout << "Bad input file: missing header: " << AdjGraphHeader << endl;
      abort();
    }

    // file consists of [type, num_vertices, num_edges, <vertex offsets>, <edges>]
    // in compressed sparse row format
    long n = atol(W[1]);
    long m = atol(W[2]);
    if (W.size() != n + m + 3) {
      cout << "Bad input file: length = "<< W.size() << " n+m+3 = " << n+m+3 << endl;
      abort(); }
    
    // tags on m at the end (so n+1 total offsets)
    sequence<intT> offsets(n+1, [&] (size_t i) {
	return (i == n) ? m : atol(W[i+3]);});
    sequence<intT> edges(m, [&] (size_t i) {
	return atol(W[n+i+3]);});

    return graph<intT>(std::move(offsets), std::move(edges), n);
  }

  pbbs::sequence<char> mmapStringFromFile(const char *filename) {
    struct stat sb;
    int fd = open(filename, O_RDONLY);
    if (fd == -1) {
      perror("open");
      exit(-1);
    }
    if (fstat(fd, &sb) == -1) {
      perror("fstat");
      exit(-1);
    }
    if (!S_ISREG (sb.st_mode)) {
      perror("not a file\n");
      exit(-1);
    }
    char *p = static_cast<char*>(mmap(0, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0));
    if (p == MAP_FAILED) {
      perror("mmap");
      exit(-1);
    }
    if (close(fd) == -1) {
      perror("close");
      exit(-1);
    }
    size_t n = sb.st_size;
    return pbbs::sequence<char>(p, n);
  }

  // template <class intT, class intE>
  // graphC<intT, intE> readGraphCFromFile(char* fname, bool mmap=false) {

  //   pbbs::sequence<char*> W;
  //   if (mmap) {
  //     cout << "mmapping file" << endl;
  //     pbbs::sequence<char> S = mmapStringFromFile(fname);
  //     // copy to new sequence
  //     pbbs::sequence<char> bytes = S;
  //     // and unmap
  //     if (munmap(S.begin(), S.size()) == -1) {
  //       perror("munmap");
  //       exit(-1);
  //     }
  //     W = stringToWords(S);
  //     cout << "mmap'd" << endl;
  //   } else {
  //     auto S = readStringFromFile(fname);
  //     W = stringToWords(S);
  //   }

  //   if (W[0] != AdjGraphHeader) {
  //     cout << "Bad input file: missing header: " << AdjGraphHeader << endl;
  //     abort();
  //   }

  //   // num vertices, num edges, edge offsets, edge pointers
  //   long len = W.size() -1;
  //   long n = atol(W[1]);
  //   long m = atol(W[2]);
  //   if (len != n + m + 2) {
  //     cout << "Bad input file: length = "<<len<< " n+m+2 = " << n+m+2 << endl;
  //     abort();
  //   }
  //   sequence<intT> offsets(n+1, [&] (size_t i) {
  // 	return (i == n) ? m : atol(W[i+3]);});
  //   sequence<intE> edges(m, [&] (size_t i) {
  // 	return atol(W[n+i+3]);});

  //   return graphC<intT,intE>(offsets,edges,n,m);
  // }

  template <class intT>
  wghGraph<intT> readWghGraphFromFile(char* fname) {
    pbbs::sequence<char> S = readStringFromFile(fname);
    pbbs::sequence<char*> W = stringToWords(S);
    if (W[0] != WghAdjGraphHeader) {
      cout << "Bad input file" << endl;
      abort();
    }

    long n = atol(W[1]);
    long m = atol(W[2]);
    if (W.size() != n + 2*m + 3) {
      cout << "Bad input file: length = "<< W.size() << " n + 2*m + 3 = " << n+2*m+3 << endl;
      abort(); }
    
    // tags on m at the end (so n+1 total offsets)
    sequence<intT> offsets(n+1, [&] (size_t i) {
	return (i == n) ? m : atol(W[i+3]);});
    sequence<intT> edges(m, [&] (size_t i) {
	return atol(W[n+i+3]);});
    sequence<intT> weights(m, [&] (size_t i) {
	return atol(W[n+i+3+m]);});

    return graph<intT>(std::move(offsets), std::move(edges), std::move(weights), n);
  }

  void errorOut(const char* s) {
    cerr << s << endl;
    throw s;
  }
  void packInt64(int64_t x, uint8_t buf[8]) {
    uint64_t xu = x;
    for (int i = 0; i < 8; ++i)
      buf[i] = (xu >> (8 * i)) & 0xff;
  }
  int64_t unpackInt64(const uint8_t buf[8]) {
    uint64_t xu = 0;
    for (int i = 0; i < 8; ++i)
      xu |= ((uint64_t)buf[i]) << (i * 8);
    return (int64_t)xu;
  }

  void writeInt(ostream& out, char buf[8], int64_t x) {
    packInt64(x, (uint8_t*)buf);
    out.write(buf, 8);
  }
  int64_t readInt(istream& in, char buf[8]) {
    in.read(buf, 8);
    return unpackInt64((uint8_t*)buf);
  }

  // template<typename intT>
  // void writeFlowGraph(ostream& out, FlowGraph<intT> g) {
  //   char buf[8];
  //   out.write("FLOWFLOW", 8);
  //   writeInt(out, buf, g.g.n);
  //   writeInt(out, buf, g.g.m);
  //   writeInt(out, buf, g.source);
  //   writeInt(out, buf, g.sink);
  //   intT o = 0;
  //   for (intT i = 0; i < g.g.n; ++i) {
  //     writeInt(out, buf, o);
  //     o += g.g.V[i].degree;
  //   }
  //   for (intT i = 0; i < g.g.n; ++i) {
  //     wghVertex<intT>& v = g.g.V[i];
  //     for (intT j = 0; j < v.degree; ++j) {
  //       writeInt(out, buf, v.Neighbors[j]);
  //       writeInt(out, buf, v.nghWeights[j]);
  //     }
  //   }
  // }
  // template<typename intT>
  // FlowGraph<intT> readFlowGraph(istream& in) {
  //   char buf[10];
  //   in.read(buf, 8);
  //   buf[8] = 0;
  //   if (strcmp(buf, "FLOWFLOW"))
  //     errorOut("Invalid flow graph input file");
  //   intT n = readInt(in, buf);
  //   intT m = readInt(in, buf);
  //   intT S = readInt(in, buf);
  //   intT T = readInt(in, buf);
  //   intT *offset = newA(intT, n);
  //   intT* adj = newA(intT, m);
  //   intT* weights = newA(intT, m);
  //   wghVertex<intT>* v = newA(wghVertex<intT>, n);
  //   for (intT i = 0; i < n; ++i) {
  //     offset[i] = readInt(in, buf);
  //     v[i].Neighbors = adj + offset[i];
  //     v[i].nghWeights = weights + offset[i];
  //     if (i > 0)
  //       v[i - 1].degree = offset[i] - offset[i - 1];
  //   }
  //   v[n - 1].degree = m - offset[n - 1];
  //   free(offset);
  //   for (intT i = 0; i < m; ++i) {
  //     adj[i] = readInt(in, buf);
  //     weights[i] = readInt(in, buf);
  //   }
  //   return FlowGraph<intT>(wghGraph<intT>(v, n, m, adj, weights), S, T);
  // }

  // const char nl = '\n';
  // template <typename intT>
  // FlowGraph<intT> writeFlowGraphDimacs(ostream& out, FlowGraph<intT> g) {
  //   out << "c DIMACS flow network description" << nl;
  //   out << "c (problem-id, nodes, arcs)" << nl;
  //   out << "p max " << g.g.n << " " << g.g.m << nl;

  //   out << "c source" << nl;
  //   out << "n " << g.source + 1 << " s" << nl;
  //   out << "c sink" << nl;
  //   out << "n " << g.sink + 1 << " t" << nl;

  //   out << "c arc description (from, to, capacity)" << nl;

  //   for (intT i = 0; i < g.g.n; ++i) {
  //     wghVertex<intT>& v = g.g.V[i];
  //     for (intT j = 0; j < v.degree; ++j) {
  //       out << "a " << i + 1 << " " << v.Neighbors[j] + 1 << " "
  //           << v.nghWeights[j] << nl;
  //     }
  //   }
  // }

  // template<typename intT>
  // struct intWghEdge {
  //   intT from, to, w;
  // };
  // int readDimacsLinePref(istream& in, const char* expected) {
  //   char type;
  //   while (in >> type) {
  //     if (type == 'c') {
  //       while (in.peek() != EOF && in.peek() != '\n')
  //         in.ignore();
  //       in >> ws;
  //       continue;
  //     } else if (!strchr(expected, type)) {
  //       errorOut((string("Unexpected DIMACS line (expected 'c' or one of '")
  // 		  + expected + "')").c_str());
  //     }
  //     return type;
  //   }
  //   return EOF;
  // }

  // template <typename intT>
  // FlowGraph<intT> readFlowGraphDimacs(istream& in) {
  //   string tmp;
  //   intT n, m;
  //   int type = readDimacsLinePref(in, "p");
  //   if (type == EOF)
  //     errorOut("Unexpected EOF while reading DIMACS file");
  //   in >> tmp >> n >> m;
  //   intWghEdge<intT>* edges = newA(intWghEdge<intT>, m);
  //   intT edgei = 0;
  //   intT* pos = newA(intT, n + 1);
  //   intT S = -1, T = -1;
  //   while (EOF != (type = readDimacsLinePref(in, "an"))) {
  //     if (type == 'n') {
  //       intT x;
  //       char st;
  //       in >> x >> st;
  //       x--;
  //       if (st == 's') S = x;
  //       else T = x;
  //     } else { // type == 'a'
  //       intT from, to, cap;
  //       in >> from >> to >> cap;
  //       from--; to--;
  //       edges[edgei] = (intWghEdge<intT>) { from, to, cap };
  //       edgei++;
  //       pos[from + 1]++;
  //     }
  //   }
  //   if (S < 0)
  //     errorOut("No source was specified in DIMACS input file");
  //   if (T < 0)
  //     errorOut("No sink was specified in DIMACS input file");
  //   if (m != edgei)
  //     errorOut("Inconsistent edge count in DIMACS input file");
  //   intT* adj = newA(intT, m);
  //   intT* weights = newA(intT, m);
  //   wghVertex<intT>* v = newA(wghVertex<intT>, n);
  //   for (intT i = 0; i < n; ++i) {
  //     pos[i + 1] += pos[i];
  //     v[i].Neighbors = adj + pos[i];
  //     v[i].nghWeights = weights + pos[i];
  //     v[i].degree = pos[i + 1] - pos[i];
  //   }
  //   for (intT i = 0; i < m; ++i) {
  //     intT& p = pos[edges[i].from];
  //     adj[p] = edges[i].to;
  //     weights[p] = edges[i].w;
  //     p++;
  //   }
  //   free(edges);
  //   free(pos);
  //   return FlowGraph<intT>(wghGraph<intT>(v, n, m, adj, weights), S, T);
  // }
};
