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
#include "../utils/types.h"
// #include "common/time_loop.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

using namespace benchIO;

// *************************************************************
//  SOME DEFINITIONS
// *************************************************************


// *************************************************************
// Parsing code (should move to common?)
// *************************************************************

// returns a pointer and a length
std::pair<char*, size_t> mmapStringFromFile(const char* filename) {
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
  if (!S_ISREG(sb.st_mode)) {
    perror("not a file\n");
    exit(-1);
  }
  char* p =
      static_cast<char*>(mmap(0, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0));
  if (p == MAP_FAILED) {
    perror("mmap");
    exit(-1);
  }
  if (close(fd) == -1) {
    perror("close");
    exit(-1);
  }
  size_t n = sb.st_size;
  return std::make_pair(p, n);
}

auto parse_fvecs(const char* filename, int maxDeg) {
  auto [fileptr, length] = mmapStringFromFile(filename);

  // Each vector is 4 + 4*d bytes.
  // * first 4 bytes encode the dimension (as an integer)
  // * next d values are floats representing vector components
  // See http://corpus-texmex.irisa.fr/ for more details.

  int d = *((int*)fileptr);

  size_t vector_size = 4 + 4*d;
  size_t num_vectors = length / vector_size;
  // std::cout << "Num vectors = " << num_vectors << std::endl;

  parlay::sequence<Tvec_point<float>> points(num_vectors);

  parlay::sequence<int> &out_nbh = *new parlay::sequence<int>(maxDeg*num_vectors);

  parlay::parallel_for(0, num_vectors*maxDeg, [&] (size_t i){out_nbh[i] = -1; });

  parlay::parallel_for(0, num_vectors, [&] (size_t i) {
    size_t offset_in_bytes = vector_size * i + 4;  // skip dimension
    float* start = (float*)(fileptr + offset_in_bytes);
    float* end = start + d;
    points[i].id = i; 
    points[i].coordinates = parlay::make_slice(start, end);
    points[i].out_nbh = parlay::make_slice(out_nbh.begin()+maxDeg*i, out_nbh.begin()+maxDeg*(i+1));
    // points[i].new_nbh = parlay::make_slice(out_nbh.begin()+maxDeg*i, out_nbh.begin()+maxDeg*(i+1));  
  });

  return points;
}

auto parse_fvecs_with_graph(const char* filename, const char* graphname){
  auto [fileptr, length] = mmapStringFromFile(filename);
  auto [graphptr, graphlength] = mmapStringFromFile(filename);

  int d = *((int*)fileptr);
  size_t vector_size = 4 + 4*d;
  size_t num_vectors = length / vector_size;

  int maxDeg = *((int*)graphptr+4);
  int num_points = *((int*)graphptr);

  if(num_vectors != num_points){
    std::cout << "ERROR: graph and data files do not match" << std::endl;
    abort();
  }

  parlay::sequence<Tvec_point<float>> points(num_vectors);

  parlay::parallel_for(0, num_vectors, [&] (size_t i) {
    points[i].id = i; 
    
    size_t offset_in_bytes = vector_size * i + 4;  // skip dimension
    float* start = (float*)(fileptr + offset_in_bytes);
    float* end = start + d;
    points[i].coordinates = parlay::make_slice(start, end);

    int* start_graph = (int*)(fileptr + 8 + maxDeg*i);
    int* end_graph = start_graph + maxDeg;
    points[i].out_nbh = parlay::make_slice(start_graph, end_graph);  
  });

  return std::make_pair(points, maxDeg);
}

auto parse_ivecs(const char* filename) {
  auto [fileptr, length] = mmapStringFromFile(filename);

  // Each vector is 4 + 4*d bytes.
  // * first 4 bytes encode the dimension (as an integer)
  // * next d values are floats representing vector components
  // See http://corpus-texmex.irisa.fr/ for more details.

  int d = *((int*)fileptr);

  size_t vector_size = 4 + 4*d;
  size_t num_vectors = length / vector_size;  

  parlay::sequence<ivec_point> points(num_vectors);

  parlay::parallel_for(0, num_vectors, [&] (size_t i) {
    size_t offset_in_bytes = vector_size * i + 4;  // skip dimension
    int* start = (int*)(fileptr + offset_in_bytes);
    int* end = start + d;
    points[i].id = i; 
    points[i].coordinates = parlay::make_slice(start, end);
  });

  return points;
}

auto parse_bvecs(const char* filename, int maxDeg) {

  auto [fileptr, length] = mmapStringFromFile(filename);
  // Each vector is 4 + d bytes.
  // * first 4 bytes encode the dimension (as an integer)
  // * next d values are unsigned chars representing vector components
  // See http://corpus-texmex.irisa.fr/ for more details.

  int d = *((int*)fileptr);
  size_t vector_size = 4 + d;
  size_t num_vectors = length / vector_size;

  parlay::sequence<Tvec_point<uint8_t>> points(num_vectors);

  parlay::sequence<int> &out_nbh = *new parlay::sequence<int>(maxDeg*num_vectors);

  parlay::parallel_for(0, num_vectors*maxDeg, [&] (size_t i){out_nbh[i] = -1;});

  parlay::parallel_for(0, num_vectors, [&] (size_t i) {
    size_t offset_in_bytes = vector_size * i + 4;  // skip dimension
    uint8_t* start = (uint8_t*)(fileptr + offset_in_bytes);
    uint8_t* end = start + d;
    points[i].id = i; 
    points[i].coordinates = parlay::make_slice(start, end);
    points[i].out_nbh = parlay::make_slice(out_nbh.begin()+maxDeg*i, out_nbh.begin()+maxDeg*(i+1));   
  });

  return points;
}

auto parse_bvecs_with_graph(const char* filename, const char* graphname) {

  auto [fileptr, length] = mmapStringFromFile(filename);
  auto [graphptr, graphlength] = mmapStringFromFile(filename);
  // Each vector is 4 + d bytes.
  // * first 4 bytes encode the dimension (as an integer)
  // * next d values are unsigned chars representing vector components
  // See http://corpus-texmex.irisa.fr/ for more details.

  int d = *((int*)fileptr);
  size_t vector_size = 4 + d;
  size_t num_vectors = length / vector_size;

  int maxDeg = *((int*)graphptr+4);
  int num_points = *((int*)graphptr);

  parlay::sequence<Tvec_point<uint8_t>> points(num_vectors);

  parlay::parallel_for(0, num_vectors, [&] (size_t i) {
    points[i].id = i; 

    size_t offset_in_bytes = vector_size * i + 4;  // skip dimension
    uint8_t* start = (uint8_t*)(fileptr + offset_in_bytes);
    uint8_t* end = start + d;
    points[i].coordinates = parlay::make_slice(start, end);

    int* start_graph = (int*)(fileptr + 8 + maxDeg*i);
    int* end_graph = start_graph + maxDeg;
    points[i].out_nbh = parlay::make_slice(start_graph, end_graph);   
  });

  return std::make_pair(maxDeg, points);
}

//graph file format begins with number of points N, then max degree
//then N+1 offsets indicating beginning and end of each vector
//then the IDs in the vector
//assumes user correctly matches data file and graph file
template<typename T>
void write_graph(parlay::sequence<Tvec_point<T>*> &v, char* outFile, int maxDeg){
  parlay::sequence<int> preamble = {static_cast<int>(v.size()), maxDeg};
  int* preamble_data = preamble.begin();
  int* graph_data = v[0]->out_nbh.begin();
  std::ofstream writer;
  writer.open(outFile, std::ios::binary | std::ios::out);
  writer.write((char *) preamble_data, 2*sizeof(int));
  writer.write((char *) graph_data, v.size()*maxDeg*sizeof(int));
  writer.close();
}


