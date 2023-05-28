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

#include <iostream>
#include "float.h"
#include <algorithm>
#include <cstring>
#include <set>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/geometry.h"
#include "common/geometryIO.h"
#include "common/parse_command_line.h"
using namespace benchIO;

using coord = double;
using point2 = point2d<coord>;
using point3 = point3d<coord>;

template<typename point>
int checkRange(parlay::sequence<int> &neighbors, parlay::sequence<point> &P,
		   double rad, int r) {
  size_t n = P.size();
  size_t num = neighbors[0];

  int wrong=0;
  int* start = neighbors.begin()+1;
  int* end = neighbors.begin()+1+num;

  parlay::slice<int*, int*> num_results = parlay::make_slice(start, end);
  auto [offsets, total] = parlay::scan(num_results);
  offsets.push_back(total);
  parlay::slice<int*, int*> ids = parlay::make_slice(end, neighbors.end());
  // std::cout << "Computed offsets" << std::endl;
  for (int j = 0; j < r; j++) {
    // std::cout << "checking point " << j << std::endl;
    int jj = parlay::hash64(j) % n;

    auto distances = parlay::tabulate(n, [&] (size_t i) {
        if(i==jj) return std::make_pair(i, DBL_MAX);
        else return std::make_pair(i, (P[jj] - P[i]).Length());
      });

    auto f = [&] (std::pair<size_t, double> a) {
      return a.second <= rad;
    };

    auto range_answers = parlay::filter(distances, f);

    std::set<int> reported_answers;
    for(int l=offsets[jj]; l<offsets[jj+1]; l++) reported_answers.insert(ids[l]);
    // std::cout << "Point " << j << " has " << reported_answers.size() << " results in range" << std::endl;
    for(auto a : range_answers){
      if(reported_answers.find(a.first) == reported_answers.end()) {wrong++; break;}
    }
    
    }
  
  if(wrong > 0) std::cout << wrong << " out of " << r << " points return wrong answer" << std::endl;
  return 0;
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,
		"[-k {1,...,100}] [-rad <rad>] [-r <numtests>] <inFile> <outfile>");
  pair<char*,char*> fnames = P.IOFileNames();
  char* iFile = fnames.first;
  char* oFile = fnames.second;

  // number of random points to test
  int r = P.getOptionIntValue("-r",10);
  double rad = P.getOptionDoubleValue("-rad", .01);
  int d = P.getOptionIntValue("-d",2);
  if (d < 2 || d > 3) P.badArgument();
  parlay::sequence<int> neighbors = readIntSeqFromFile<int>(oFile);
  if (d == 2) {
    parlay::sequence<point2> PIn = readPointsFromFile<point2>(iFile);
    return checkRange(neighbors, PIn, rad, r);
  } else if (d==3){
    parlay::sequence<point3> PIn = readPointsFromFile<point3>(iFile);
    return checkRange(neighbors, PIn, rad, r);
  }

}
