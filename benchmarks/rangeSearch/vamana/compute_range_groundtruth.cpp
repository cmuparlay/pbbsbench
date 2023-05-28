#include <iostream>
#include <algorithm>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/io.h"
#include "../utils/types.h"
#include "../utils/NSGDist.h"
#include "../utils/parse_files.h"
// #include "../bench/benchUtils.h"

using pid = std::pair<int, float>;


template<typename T>
void compute_groundtruth(parlay::sequence<Tvec_point<T>> &B, parlay::sequence<Tvec_point<T>> &Q, double radius,
    parlay::sequence<parlay::sequence<pid>> &result){
    
    unsigned d = (B[0].coordinates).size();
    parlay::parallel_for(0, Q.size(), [&] (size_t i){
        parlay::sequence<pid> in_range;
        for(size_t j=0; j<B.size(); j++){
            float dist = distance((Q[i].coordinates).begin(), (B[j].coordinates).begin(), d);
            if(dist < radius) in_range.push_back(std::make_pair(j, dist));
        }
        result[i] = in_range;
    });
}

void write_rangeres(parlay::sequence<parlay::sequence<pid>> &result, const std::string outFile){
    // std::cout << "Writing file with dimension " << result[0].size() << std::endl;
    std::cout << "File contains groundtruth for " << result.size() << " data points" << std::endl;

    // auto less = [&] (pid a, pid b) {return a.second < b.second;};

    int n = result.size();

    int num_matches = 0;
    for(int i=0; i<result.size(); i++) num_matches += result[i].size();

    parlay::sequence<int> preamble = {n, num_matches};

    auto offsets = parlay::tabulate(result.size(), [&] (size_t i){
      return static_cast<int>(result[i].size());
    });
  
    auto ids = parlay::tabulate(result.size(), [&] (size_t i){
        parlay::sequence<int> data;
        if(result[i].size()==0) return data;
        else{
          for(int j=0; j<result[i].size(); j++){data.push_back(result[i][j].first);}
          return data;
        }  
    });

    auto ids_data = parlay::flatten(ids);

    auto distances = parlay::tabulate(result.size(), [&] (size_t i){
        parlay::sequence<int> data;
        if(result[i].size()==0) return data;
        else{
          for(int j=0; j<result[i].size(); j++){data.push_back(static_cast<int>(result[i][j].second));}
          return data;
        }       
    });

    auto distances_data = parlay::flatten(distances);

    int* preamble_data = preamble.begin();
    int* offset_data = offsets.begin();
    int* id_data = ids_data.begin();
    int* distance_data = distances_data.begin();

    std::ofstream writer;
    writer.open(outFile, std::ios::binary | std::ios::out);
    writer.write((char *) preamble_data, 2*sizeof(int));
    writer.write((char *) offset_data, result.size()*sizeof(int));
    writer.write((char *) id_data, ids_data.size()*sizeof(int));
    writer.write((char *) distance_data, distances_data.size()*sizeof(int));
    writer.close();
}


int main(int argc, char* argv[]) {
  if (argc != 5) {
    std::cout << "usage: compute_range_groundtruth <base> <query> <radius> <oFile>" << std::endl;
    return 1;
  }
  double rad = std::atof(argv[3]);
  std::cout << "Computing all points within radius " << rad << " of each query point" << std::endl;
  bool fvecs = true;
  std::string filename = std::string(argv[1]);
  std::string::size_type n = filename.size();
  if(filename[n-5] == 'b') fvecs = false;

  int maxDeg = 0;
  auto [md, B] = parse_uint8bin(argv[1], NULL, maxDeg);
  auto [fd, Q] = parse_uint8bin(argv[2], NULL, maxDeg);
  parlay::sequence<parlay::sequence<pid>> result(Q.size());
  compute_groundtruth(B, Q, rad, result);


  write_rangeres(result, std::string(argv[4]));

  return 0;
}