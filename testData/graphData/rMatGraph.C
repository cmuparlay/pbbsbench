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

#include "common/IO.h"
#include "common/graph.h"
#include "common/graphIO.h"
#include "common/graphUtils.h"
#include "common/parse_command_line.h"
#include "parlay/parallel.h"
using namespace benchIO;
using namespace dataGen;
using namespace std;

template <class intV>
struct rMat {
  using edgeT = edge<intV>;
  double a, ab, abc;
  size_t n; 
  size_t h;
  rMat(size_t _n, size_t _seed, 
       double _a, double _b, double _c) {
    n = _n; a = _a; ab = _a + _b; abc = _a+_b+_c;
    h = dataGen::hash<size_t>(_seed);
    if (!(abc <= 1.0)) {
      std::cout << "in rMat: a + b + c add to more than 1" << std::endl;
      abort();
    }
    if (!((1 << parlay::log2_up(n)) == n)) {
      std::cout << "in rMat: n not a power of 2" << std::endl;
      abort();
    }		 
  }

  edgeT rMatRec(size_t nn, size_t randStart, size_t randStride) {
    if (nn==1) return edgeT(0,0);
    else {
      edgeT x = rMatRec(nn/2, randStart + randStride, randStride);
      double r = dataGen::hash<double>(randStart);
      if (r < a) return x;
      else if (r < ab) return edgeT(x.u,x.v+nn/2);
      else if (r < abc) return edgeT(x.u+nn/2, x.v);
      else return edgeT(x.u+nn/2, x.v+nn/2);
    }
  }

  edge<intV> operator() (size_t i) {
    size_t randStart = dataGen::hash<size_t>((2*i)*h);
    size_t randStride = dataGen::hash<size_t>((2*i+1)*h);
    return rMatRec(n, randStart, randStride);
  }
};

template <class intV>
edgeArray<intV> edgeRmat(size_t n, size_t m, size_t seed, 
			 double a, double b, double c) {
  size_t nn = (1 << parlay::log2_up(n));
  rMat<intV> g(nn,seed,a,b,c);
  auto E = parlay::tabulate(m, [&] (size_t i) -> edge<intV> {return g(i);});
  return edgeArray<intV>(std::move(E), nn, nn);
}


int main(int argc, char* argv[]) {
  commandLine P(argc,argv,
		"[-m <numedges>] [-s <intseed>] [-o] [-j] [-a <a>] [-b <b>] [-c <c>] n <outFile>");
  pair<size_t,char*> in = P.sizeAndFileName();
  size_t n = in.first;
  char* fname = in.second;
  double a = P.getOptionDoubleValue("-a",.5);
  double b = P.getOptionDoubleValue("-b",.1);
  double c = P.getOptionDoubleValue("-c", b);
  size_t m = P.getOptionLongValue("-m", 10*n);
  size_t seed = P.getOptionLongValue("-s", 1);
  bool adjArray = P.getOption("-j");
  bool ordered = P.getOption("-o");

  edgeArray<size_t> EA = edgeRmat<size_t>(n, m, seed, a, b, c);
  writeGraphFromEdges(EA, fname, adjArray, ordered);
}
