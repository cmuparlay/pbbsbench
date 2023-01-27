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
#include "../../ANN/bench/benchUtils.h"
// #include "common/time_loop.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

using namespace benchIO;

auto parse_uint8bin(const char* filename, int maxDeg){
    auto [fileptr, length] = mmapStringFromFile(filename);
    int num_vectors = *((int*) fileptr);
    int d = *((int*) (fileptr+4));
    std::cout << "Detected " << num_vectors << " points with dimension " << d << std::endl;
    parlay::sequence<Tvec_point<uint8_t>> points(num_vectors);

    parlay::sequence<int> &out_nbh = *new parlay::sequence<int>(maxDeg*num_vectors);

    parlay::parallel_for(0, num_vectors*maxDeg, [&] (size_t i){out_nbh[i] = -1; });

    parlay::parallel_for(0, num_vectors, [&] (size_t i) {
        uint8_t* start = (uint8_t*)(fileptr + 8 + i*d); //8 bytes at the start for size + dimension
        uint8_t* end = start + d;
        points[i].id = i; 
        points[i].coordinates = parlay::make_slice(start, end);
        points[i].out_nbh = parlay::make_slice(out_nbh.begin()+maxDeg*i, out_nbh.begin()+maxDeg*(i+1));
    });

    return points;
}

auto parse_uint8bin_with_graph(const char* filename, const char* graphname){
    auto [fileptr, length] = mmapStringFromFile(filename);
    auto [graphptr, graphlength] = mmapStringFromFile(graphname);

    int num_vectors = *((int*) fileptr);
    int d = *((int*) (fileptr+4));

    int maxDeg = *((int*)(graphptr+4));
    int num_points = *((int*)graphptr);

    if(num_vectors != num_points){
        std::cout << "ERROR: graph and data files do not match" << std::endl;
        abort();
    }

    std::cout << "Detected " << num_vectors << " points with dimension " << d << std::endl;
    parlay::sequence<Tvec_point<uint8_t>> points(num_vectors);

    parlay::parallel_for(0, num_vectors, [&] (size_t i) {
        points[i].id = i; 

        uint8_t* start = (uint8_t*)(fileptr + 8 + i*d); //8 bytes at the start for size + dimension
        uint8_t* end = start + d;
        points[i].coordinates = parlay::make_slice(start, end);

        int* start_graph = (int*)(graphptr + 8 + 4*maxDeg*i);
        int* end_graph = start_graph + maxDeg;
        points[i].out_nbh = parlay::make_slice(start_graph, end_graph);  
    });

    return std::make_pair(maxDeg, points);
}

auto parse_rangeres(const char* filename){
    auto [fileptr, length] = mmapStringFromFile(filename);
    int num_points = *((int*) fileptr);
    int num_matches = *((int*) (fileptr+4));
    
    std::cout << "Detected " << num_points << " query points with " << num_matches << " unique matches" << std::endl;
    int* start = (int*)(fileptr+8);
    int* end = start + num_points;
    parlay::slice<int*, int*> num_results = parlay::make_slice(start, end);
    auto [offsets, total] = parlay::scan(num_results);
    offsets.push_back(total);
    parlay::sequence<ivec_point> points(num_points);

    auto id_offset = 4*num_points+8;
    auto dist_offset = id_offset + 4*num_matches; 
    parlay::parallel_for(0, num_points, [&] (size_t i) {
        int* start = (int*)(fileptr + id_offset + 4*offsets[i]); 
        int* end = (int*)(fileptr + id_offset + 4*offsets[i+1]);
        // float* dist_start = (float*)(fileptr + dist_offset + 4*offsets[i]); 
        // float* dist_end = (float*)(fileptr + dist_offset + 4*offsets[i+1]);
        points[i].id = i;
        points[i].coordinates = parlay::make_slice(start, end);
        // points[i].distances = parlay::make_slice(dist_start, dist_end);
    });

    // int i=0;
    // int j=0;
    // int cutoff = 10;
    // while(j<cutoff){
    //     if(points[i].coordinates.size() > 0){
    //         std::cout << "The " << i << "th query has " << (points[i].coordinates).size() << " matches" << std::endl;
    //         std::cout << "The match IDs are " ;
    //         for(int k=0; k<points[i].coordinates.size(); k++){
    //             std::cout << points[i].coordinates[k] << ", " ;
    //         }
    //         std::cout << std::endl;
    //         std::cout << "The match distances are " ;
    //         for(int k=0; k<points[i].distances.size(); k++){
    //             std::cout << points[i].distances[k] << ", " ;
    //         }
    //         std::cout << std::endl;
    //         j++;
    //     }
    //     i++;
    // }
    // for(int i=0; i<10; i++){
    //     std::cout << "The " << i << "th query has " << (points[i].coordinates).size() << " matches" << std::endl;
    // }

    return points;
}