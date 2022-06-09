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
#include <algorithm>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/geometry.h"
#include "common/geometryIO.h"
#include "common/parse_command_line.h"
#include "common/time_loop.h"
#include "benchUtils.h"


#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

using namespace benchIO;

bool report_stats = true;


// *************************************************************
//  TIMING
// *************************************************************

template<typename T>
void timeNeighbors(parlay::sequence<Tvec_point<T>> &pts,
  int rounds, int maxDeg, int beamSize, double delta, double alpha, char* outFile) 
{
  size_t n = pts.size();
  auto v = parlay::tabulate(n, [&] (size_t i) -> Tvec_point<T>* {
      return &pts[i];});
 
  time_loop(rounds, 0, 
  [&] () {},
  [&] () {
    ANN<T>(v, maxDeg, beamSize, alpha, delta);
  },
  [&] () {});  

}

template<typename T>
void timeNeighbors(parlay::sequence<Tvec_point<T>> &pts, parlay::sequence<Tvec_point<T>> &qpoints,
  int k, int rounds, int maxDeg, int beamSize, int beamSizeQ, double delta, double alpha, char* outFile) 
{
  size_t n = pts.size();
  auto v = parlay::tabulate(n, [&] (size_t i) -> Tvec_point<T>* {
      return &pts[i];});

  size_t q = qpoints.size();
  auto qpts =  parlay::tabulate(q, [&] (size_t i) -> Tvec_point<T>* {   
      return &qpoints[i];});

    time_loop(rounds, 0, 
      [&] () {},
      [&] () {
        ANN<T>(v, k, maxDeg, beamSize, beamSizeQ, alpha, delta, qpts); 
      },
      [&] () {});

    if (outFile != NULL) {
      int m = q * (k+1);
      parlay::sequence<int> Pout(m);
      parlay::parallel_for (0, q, [&] (size_t i) {
        Pout[(k+1)*i] = qpts[i]->id;
        for (int j=0; j < k; j++)
          Pout[(k+1)*i + j+1] = (qpts[i]->ngh)[j];
      });
      writeIntSeqToFile(Pout, outFile);
    }
  

}

// Infile is a file in .fvecs format
int main(int argc, char* argv[]) {
    commandLine P(argc,argv,
    "[-a <alpha>] [-d <delta>] [-R <deg>]"
        "[-L <bm>] [-k <k> ] [-Q <bmq>] [-q <qF>] [-o <oF>] [-r <rnds>] [-b <algoOpt>] <inFile>");

  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  char* qFile = P.getOptionValue("-q");
  int R = P.getOptionIntValue("-R", 5);
  if (R < 1) P.badArgument();
  int L = P.getOptionIntValue("-L", 10);
  if (L < 1) P.badArgument();
  int Q = P.getOptionIntValue("-Q", L);
  if (Q < 1) P.badArgument();
  int rounds = P.getOptionIntValue("-r", 1);
  int k = P.getOptionIntValue("-k", 1);
  if (k > 1000 || k < 1) P.badArgument();
  double alpha = P.getOptionDoubleValue("-a", 1.2);
  double delta = P.getOptionDoubleValue("-d", .01);
  int HCNNG = P.getOptionIntValue("-b", 0);

  bool fvecs = true;
  std::string filename = std::string(iFile);
  std::string::size_type n = filename.size();
  if(filename[n-5] == 'b') fvecs = false;

  int maxDeg;
  if(HCNNG != 0) maxDeg = R*L;
  else maxDeg = R;



  if(fvecs){ //vectors are floating point coordinates
    parlay::sequence<Tvec_point<float>> points = parse_fvecs(iFile, maxDeg);
    if(qFile != NULL){
      parlay::sequence<Tvec_point<float>> qpoints = parse_fvecs(qFile, 0);
      timeNeighbors<float>(points, qpoints, k, rounds, R, L, Q, delta, alpha, oFile);
    }
    else timeNeighbors<float>(points, rounds, R, L, delta, alpha, oFile);
  } 
  else{ //vectors are uint8 coordinates
    parlay::sequence<Tvec_point<uint8_t>> points = parse_bvecs(iFile, maxDeg);
    if(qFile != NULL){
      parlay::sequence<Tvec_point<uint8_t>> qpoints = parse_bvecs(qFile, 0);
      timeNeighbors<uint8_t>(points, qpoints, k, rounds, R, L, Q, delta, alpha, oFile);
    }
    else timeNeighbors<uint8_t>(points, rounds, R, L, delta, alpha, oFile);}  
}



//REGULAR CORRECTNESS CHECK
// if (outFile != NULL) {
      // int m = q * (k+1);
      // parlay::sequence<int> Pout(m);
      // parlay::parallel_for (0, q, [&] (size_t i) {
      //   Pout[(k+1)*i] = qpts[i]->id;
      //   for (int j=0; j < k; j++)
      //     Pout[(k+1)*i + j+1] = (qpts[i]->ngh)[j];
      // });
      // writeIntSeqToFile(Pout, outFile);
    // }


//DIRECTED DEGREES
// if (outFile != NULL) {
      // parlay::sequence<int> degrees(n);
      // parlay::parallel_for(0, n, [&] (size_t i){
      //   degrees[i] = v[i]->out_nbh.size();
      // });
      // writeIntSeqToFile(degrees, outFile);
    // }

//COORDINATES OF FILE
// if (outFile != NULL) {
      // int d = v[0]->coordinates.size();
      // int m = n*d; //total number of ints in the file
      // parlay::sequence<int> vals(m);
      // parlay::parallel_for(0, n, [&] (size_t i){
          // for(int j=0; j<d; j++){
          //   vals[i*d+j] = (int) v[i]->coordinates[j];
          // }
      // });
      // writeIntSeqToFile(vals, outFile);
    // }

//UNDIRECTED DEGREES (REALLY SLOW AND BAD)
 // if (outFile != NULL) {
 //      parlay::sequence<int> udegrees(n);
 //      parlay::sequence<parlay::sequence<index_pair>> to_flatten = parlay::sequence<parlay::sequence<index_pair>>(n);
 //      parlay::parallel_for(0, n, [&] (size_t i){
 //        size_t m = v[i]->out_nbh.size();
 //        parlay::sequence<index_pair> edges = parlay::sequence<index_pair>(2*m);
 //        for(int j=0; j<m; j++){
 //          edges[2*j] = std::make_pair(i, v[i]->out_nbh[j]);
 //          edges[2*j+1] = std::make_pair(v[i]->out_nbh[j], i);
 //        }
 //        to_flatten[i] = edges;
 //      });
 //      std::cout << "here1" << std::endl; 
 //    auto edges_unsorted = parlay::flatten(to_flatten);
 //    std::cout << "here2" << std::endl; 
 //    auto grouped_edges = parlay::group_by_key(edges_unsorted);
 //    std::cout << "here3" << std::endl; 
 //    parlay::parallel_for(0, n, [&] (size_t i){
 //      int count = 0;
 //      // parlay::sequence<int> edge_ids = grouped_edges[i].second;
 //      std::cout << grouped_edges[i].second.size() << std::endl; 
 //      for(int j=0; j<grouped_edges[i].second.size(); j++){
 //        count+=1;
        
 //        if(j<grouped_edges[i].second.size()-1){
 //          int current = grouped_edges[i].second[j];
 //          int next = grouped_edges[i].second[j+1];
 //          if(current == next) j+=1;
 //        }
        
 //      }
 //      udegrees[i] = count;
 //    });
 //    auto sortedDegrees = parlay::sort(udegrees);
 //    writeIntSeqToFile(sortedDegrees, outFile);
 //    }

  //GRAPH FORMAT

   // if (outFile != NULL) {
   //    parlay::sequence<int> graph(n*(maxDeg+1));
   //    parlay::parallel_for(0, n, [&] (size_t i){
   //      graph[i*(maxDeg+1)] = (int) i; 
   //      int degree = v[i]->out_nbh.size();
   //      for(int j=0; j<degree; j++) graph[i*(maxDeg+1)+1+j] = v[i]->out_nbh[j]; 
   //      int rem = maxDeg - degree;
   //      for(int j=0; j<rem; j++) graph[i*(maxDeg+1)+1+degree+j] = -1;
   //    });
   //    writeIntSeqToFile(graph, outFile);
   //  }




