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
int algorithm_version = 0;
// 0=root based, 1=bit based

#include <algorithm>
#include <math.h> 
#include <queue>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/geometry.h"
#include "k_nearest_neighbors.h"

// find the k nearest neighbors for all points in tree
// places pointers to them in the .ngh field of each vertex
template <int max_k, class vtx>
void ANN(parlay::sequence<vtx*> &v, int k) {
  timer t("ANN",report_stats);

  {
    using knn_tree = k_nearest_neighbors<vtx, max_k>;
    using point = typename knn_tree::point;
    using node = typename knn_tree::node;
    using box = typename knn_tree::box;
    using box_delta = std::pair<box, double>;
  
    //calculate bounding box around the whole point set
    box whole_box = knn_tree::o_tree::get_box(v);

    //split initial vertices into two sequences: one to build the tree with
    //and one to later insert point by point
    size_t n = v.size(); 
    size_t init = 3*n/4;
    size_t ins = n-init;
    parlay::sequence<vtx*> v1(init);
    parlay::sequence<vtx*> v2(ins);
    parlay::parallel_for(0, n, [&] (size_t i){
      if(i<init) v1[i] = v[i];
      else v2[i-init] = v[i];
    }, 1
    );

    //build tree with bounding box
    knn_tree T(v1, whole_box);
    t.next("build tree");

    //prelims for insert  
    int dims = v[0]->pt.dimension();
    box_delta bd = T.get_box_delta(dims);


    // for(int j = 0; j < v2.size(); j++)
    //   T.insert_point(v2[j], bd.first, bd.second); 
    
    parlay::parallel_for(0, parlay::num_workers(), [&] (size_t i) {
      for(int j = i; j < v2.size(); j+=parlay::num_workers()) {
        T.insert_point(v2[j], bd.first, bd.second); 
      }
      // std::cout << "thread finished inserting\n";
    }, 1, true);

    t.next("insert points");

    if (report_stats) 
      std::cout << "depth = " << T.tree.load()->depth() << std::endl;

    if (algorithm_version == 0) { // this is for starting from root 

      // find nearest k neighbors for each point
      size_t n = v.size();
      parlay::parallel_for (0, n, [&] (size_t i) {
	       T.k_nearest(v[i], k);
      }, 1);
    
    } else {

        int dims = (v[0]->pt).dimension();  
        node* root = T.tree.load(); 
        box_delta bd = T.get_box_delta(dims);
        size_t n = v.size();
        parlay::parallel_for(0, n, [&] (size_t i) {
          T.k_nearest_leaf(v[i], T.find_leaf(v[i]->pt, root, bd.first, bd.second), k);
        }
        );
    }

    t.next("try all");
    if (report_stats) {
      auto s = parlay::delayed_seq<size_t>(v.size(), [&] (size_t i) {return v[i]->counter;});
      size_t i = parlay::max_element(s) - s.begin();
      size_t sum = parlay::reduce(s);
      std::cout << "max internal = " << s[i] 
		<< ", average internal = " << sum/((double) v.size()) << std::endl;
      t.next("stats");
    }
    t.next("delete tree");   


};
}

