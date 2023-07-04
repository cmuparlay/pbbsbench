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

bool report_stats = true;

#include <algorithm>
#include <math.h> 
#include <queue>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/geometry.h"
#include "../../concurrentKNN/octTree/k_nearest_neighbors.h"
#include "../../concurrentKNN/octTree/rand_r_32.h"

/* compile with debug flags:
g++-11 -DHOMEGROWN -pthread -mcx16 -g -std=c++17 -I .  -include neighbors.h -o neighbors_debug ../bench/neighborsTime.C -fsanitize=address -fsanitize=undefined -DHOMEGROWN -pthread -ldl -L/usr/local/lib -ljemalloc

ASAN_OPTIONS=detect_leaks=0 PARLAY_NUM_THREADS=72 numactl -i all ./neighbors_debug -d 2 -k 1 -o oFile ../geometryData/data/2DinCube_10M
*/

// find the k nearest neighbors for all points in tree
// places pointers to them in the .ngh field of each vertex
template <class vtx>
void RANGE(parlay::sequence<vtx*> &v, double rad, int p, double trial_time, int update_percent) {
  timer t("ANN",report_stats);

  {
    using knn_tree = k_nearest_neighbors<vtx, 1>;
    using point = typename knn_tree::point;
    using node = typename knn_tree::node;
    using box = typename knn_tree::box;
    using box_delta = std::pair<box, double>;
  
#ifdef Versioned
    std::cout << "using multiversioning" << std::endl;
#else
    std::cout << "without multiversioning" << std::endl;
#endif

#ifndef NoHelp
    std::cout << "using lock-free locks" << std::endl;
#else
    std::cout << "using blocking locks" << std::endl;
#endif

#ifdef HandOverHand
    std::cout << "using hand-over-hand locking" << std::endl;
#else
    std::cout << "using path locking" << std::endl;
#endif

    std::cout << "threads: " << p << std::endl;
    std::cout << "update_percent: " << update_percent << std::endl;
    std::cout << "query radius: " << rad << std::endl;
    std::cout << "trial_time: " << trial_time << std::endl;

    //calculate bounding box around the whole point set
    box whole_box = knn_tree::o_tree::get_box(v);     

    //split initial vertices into two sequences: one to build the tree with
    //and one to later insert point by point
    size_t n = v.size();
    v = parlay::random_shuffle(v);
    node_allocator<vtx>.shuffle(n); 
    size_t init = n/2;
    // size_t ins = n-init;
    parlay::sequence<vtx*> v_init(init);  
    // parlay::sequence<vtx*> v2(ins);
    parlay::parallel_for(0, init, [&] (size_t i){
      v_init[i] = v[i];
    }, 1);

    parlay::sequence<size_t> totals(p);
    parlay::sequence<long> addeds(p);
    parlay::sequence<long> ins_fails(p);
    parlay::sequence<long> del_fails(p);
    
    // t.next("setup benchmark");

    //build tree with bounding box
    knn_tree T(v_init, whole_box);
    // t.next("build tree");
   
    // run benchmark
    t.start();
    auto start = std::chrono::system_clock::now();
    std::atomic<bool> finish = false;

    parlay::sequence<parlay::sequence<int>> visited_stats(p);
    parlay::sequence<parlay::sequence<int>> size_stats(p);

    parlay::parallel_for(0, p, [&] (size_t i) {
      int cnt = 0;
      size_t total = 0;
      long added = 0;
      long ins_failed = 0;
      long del_failed = 0;
      my_rand::init(i);
      while (true) {
        // every once in a while check if time is over
        if (cnt == 100) {
          cnt = 0;
          auto current = std::chrono::system_clock::now();
          double duration = std::chrono::duration_cast<std::chrono::seconds>(current - start).count();
          if (duration > trial_time || finish) {
            totals[i] = total;
            addeds[i] = added;
            ins_fails[i] = ins_failed;
            del_fails[i] = del_failed;
            return;
          }
        }
        int op_type = my_rand::get_rand()%100;
        int idx = my_rand::get_rand()%n;
        // std::cout << op_type << " " << idx << std::endl;
        if (op_type < update_percent/2) { 
          if(T.insert_point(v[idx])) added++;
          else ins_failed++; 
        } else if (op_type < update_percent) { 
          if(T.delete_point(v[idx])) added--;
          else del_failed++; 
        } else { 
            auto ans = T.range_search(v[idx], rad);
            visited_stats[i].push_back(std::get<1>(ans)); 
            size_stats[i].push_back(static_cast<int>(std::get<2>(ans).size()));
        }
        cnt++;
        total++;
      }
    }, 1);
    double duration = t.stop();

    //std::cout << duration << " : " << trial_time << std::endl;
    size_t num_ops = parlay::reduce(totals);
    std::cout << "throughput (Mop/s): "
            << num_ops / (duration * 1e6) << std::endl << std::endl;

    parlay::sequence<int> s = parlay::flatten(visited_stats);
    size_t i = parlay::max_element(s) - s.begin();
    size_t sum = parlay::reduce(s);
    std::cout << "max internal = " << s[i] << ", average internal = " << sum/((double) s.size()) << std::endl;
    parlay::sequence<int> sizes = parlay::flatten(size_stats);
    size_t j = parlay::max_element(sizes) - sizes.begin();
    size_t sum_size = parlay::reduce(sizes);
    std::cout << "max results = " << sizes[j] << ", average results = " << sum_size/((double) sizes.size()) << std::endl;
    std::cout << "failed ins: " << parlay::reduce(ins_fails) << std::endl;
    std::cout << "failed del: " << parlay::reduce(del_fails) << std::endl;
    std::cout << "total ops: " << num_ops << endl;

    if (report_stats) {
      std::cout << "depth = " << T.tree.load()->depth() << std::endl;
    }
};
}

