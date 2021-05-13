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

// Code roughly based on C4.5 decision trees
// Works with discrete or continuous features.
// Tries all binary cuts for continuous variables, and picks best.
// The target label must be discrete (i.e. no regression).
// Adds cost of encoding a node to each node for the purpose of pruning.

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <set>
#include <iterator>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/io.h"
#include "parlay/internal/get_time.h"
#include "classify.h"
#include <random>

using namespace parlay;
double infinity = std::numeric_limits<double>::infinity();

// some parameters

// minimum size of a node
size_t min_size = 1;
// the following helps prevent overfitting if > 0.0
// probably should be in range [0..1]
double encode_node_factor = 0.0; 

// random forest parameters
const bool use_random_forest = false;
const bool use_random_entries_sample = false;
const double random_entries_sample_factor = 1.0;
const bool use_random_features_sample = false;
const double random_features_sample_factor = 0.5;
const size_t num_trees = 11;

struct tree {
  bool is_leaf;
  int feature_index;
  int feature_cut;
  int best;
  size_t size;
  sequence<tree*> children;
  tree(int i, int c, int best, sequence<tree*> children)
    : is_leaf(false), best(best), feature_index(i), feature_cut(c), children(children) {
    size = reduce(map(children, [] (tree* t) {return t->size;}));}
  tree(int best) : is_leaf(true), best(best), size(1) {}
};

tree* Leaf(int best) {
  if (best > max_value) abort();
  return new tree(best);
}

tree* Internal(int i, int cut, int majority, sequence<tree*> children) {
  return new tree(i, cut, majority, children);
}

// To be put into parlay

auto histogramSequential(auto &s, int max){
  sequence<int> freqs = tabulate(max, [&] (size_t i) { return 0; });
  for(int i = 0; i < s.size(); i++) {
    freqs[s[i]]++;
  }
  return freqs; 
}

template <typename S1, typename S2>
auto delayed_zip(S1 const &a, S2 const &b) {
  return delayed_tabulate(a.size(), [&] (size_t i) {return std::pair(a[i],b[i]);});
}

template <typename S1>
auto all_equal(S1 const &a) {
  return (a.size() == 0) || (count(a, a[0]) == a.size());
}

template <typename S1>
auto majority(S1 const &a, size_t m) { 
  auto x = histogramSequential(a,m);
  return max_element(x) - x.begin();
}

// End to be put into parlay

// entropy * total given a histogram a (total = sum(a))
template <typename Seq>
double entropy(Seq a, int total) {
  double ecost = encode_node_factor * log2(float(1 + total)); // to prevent overfitting
  return ecost + reduce(map(a, [=] (int l) {
      return (l > 0) ? -(l * log2(float(l)/total)) : 0.0;}));
}

auto cond_info_continuous(feature const &a, feature const &b) {
  int num_buckets = a.num * b.num;
  size_t n = a.vals.size();
  auto sums = histogram(delayed_tabulate(n, [&] (size_t i) {
			   return a.vals[i] + b.vals[i]*a.num;}), num_buckets);
  sequence<int> low_counts(a.num, 0);
  sequence<int> high_counts(a.num, 0);
  for (int i=0; i < b.num; i++) 
    for (int j=0; j < a.num; j++) high_counts[j] += sums[a.num*i + j];
  double cur_e = infinity;
  int cur_i = 0;
  int m = 0;
  for (int i=0; i < b.num-1; i++) {
    for (int j=0; j < a.num; j++) {
      low_counts[j] += sums[a.num*i + j];
      high_counts[j] -= sums[a.num*i + j];
      m += sums[a.num*i + j];
    }
    double e = entropy(low_counts, m) + entropy(high_counts, n - m);
    if (e < cur_e) {
      cur_e = e;
      cur_i = i+1;
    }
  }
  return std::pair(cur_e, cur_i);
}

// information content of s (i.e. entropy * size)
double info(row s, int num_vals) {
  size_t n = s.size();
  if (n == 0) return 0.0;
  auto x = histogramSequential(s, num_vals);
  return entropy(x, n);
}

// info of a conditioned on b
double cond_info_discrete(feature const &a, feature const &b) {
  int num_buckets = a.num * b.num;
  size_t n = a.vals.size();
  auto sums = histogram(delayed_tabulate(n, [&] (size_t i) {
			   return a.vals[i] + b.vals[i]*a.num;}), num_buckets);
  return reduce(tabulate(b.num, [&] (size_t i) {
      auto x = sums.cut(i*a.num,(i+1)*a.num);				      
      return entropy(x, reduce(x));}));
}

// information needed to describe the node
double node_cost(int n, int num_features, int num_groups) {
  return log2(float(num_features));
}

auto build_tree(features &A, bool verbose, auto usedFeatures) {
  int num_features = A.size();
  int num_entries = A[0].vals.size();

  int majority_value = (num_entries == 0) ? -1 : majority(A[0].vals, A[0].num);
  if (num_entries < 2 || all_equal(A[0].vals))
    return Leaf(majority_value);
  double label_info = info(A[0].vals,A[0].num);

  sequence<int> featuresIdxs = tabulate(num_features - 1, [&] (int i) { return i; }, 1);
  
  if(use_random_features_sample) {
    sequence<int> nonUsed = filter(featuresIdxs, [&] (int f) {
      return usedFeatures.find(f + 1) == usedFeatures.end();
    });
    int n = nonUsed.size();
    sequence<int> permuted = tabulate(n, [&] (int i) { return nonUsed[i]; }, 1);
    thread_local std::random_device rng;
    thread_local std::default_random_engine e(rng());
    for(int i = n - 1; i >= 1; i--) {
      std::uniform_int_distribution<int> uniform_dist(0, i);
      int j = uniform_dist(e) % (i + 1);
      std::swap(permuted[i], permuted[j]);
    }
    int num = n * random_features_sample_factor;
    if(num <= 0) num = 1;
    else if(num > n) num = n;
    featuresIdxs = tabulate(num, [&] (int i) { return permuted[i]; });
  }

  auto costs = tabulate(featuresIdxs.size(), [&] (int j) {
      int i = featuresIdxs[j];
      if (A[i+1].discrete) {
	return std::tuple(cond_info_discrete(A[0], A[i+1]), i+1, -1);
      } else {
	auto [info, cut] = cond_info_continuous(A[0], A[i+1]);
	return std::tuple(info, i+1, cut);
      }},1);
  auto min1 = [&] (auto a, auto b) {return (std::get<0>(a) < std::get<0>(b)) ? a : b;};
  auto min_m = make_monoid(min1, std::tuple(infinity, 0, 0));
  auto [best_info, best_i, cut] = reduce(costs, min_m);
  double threshold = log2(float(featuresIdxs.size()));  

  if (verbose)
    cout << num_entries << ", " << num_features << ", " << best_i << ", " << cut << ", " << label_info << ", " 
	 << best_info << ", " << featuresIdxs.size() << ", " << threshold << endl;

  if (label_info - best_info < threshold){
    return Leaf(majority_value);
  } else {
    int m;
    row split_on;
    if (A[best_i].discrete) {
      m = A[best_i].num;
      split_on = A[best_i].vals;
    } else {
      m = 2;
      split_on =  map(A[best_i].vals, [&] (value x) -> value {return x >= cut;});
    }

    features F = map(A, [&] (feature a) {return feature(a.discrete, a.num);});
    sequence<features> B(m, F);
    parallel_for (0, num_features, [&] (size_t j) {
      auto x = group_by_index(delayed_zip(split_on, A[j].vals), m);
      for (int i=0; i < m; i++) B[i][j].vals = std::move(x[i]);
    }, 1);
    //A.clear();

    if(use_random_features_sample && A[best_i].discrete)
    	usedFeatures.insert(best_i);
    auto children = map(B, [&] (features &a) {return build_tree(a, verbose, usedFeatures);}, 1);

    return Internal(best_i - 1, cut, majority_value, children); //-1 since first is label
  }
}

int classify_row(tree* T, row const&r) {
  if (T->is_leaf) {
    return T->best;
  } else if (T->feature_cut == -1) { // discrete partition
    // could be a feature value in the test data that did not appear in training data
    // in this case return the best
    if (!(r[T->feature_index] < T->children.size())) return T->best; // -1;
    // go to child based on feature value
    int val = classify_row(T->children[r[T->feature_index]], r);
    return (val == -1) ? T->best : val;
  } else {  // continuous cut
    // go to child based on whether below or at-above cut value
    int idx = (r[T->feature_index] < T->feature_cut) ? 0 : 1;
    int val = classify_row(T->children[idx], r);
    return (val == -1) ? T->best : val;
  }
}

tree* buildRandomTree(features const &Train, bool verbose) {
  features A = Train;
 
  int num_features = A.size();
  int num_entries = A[0].vals.size();

  features B = Train;  
  if(use_random_entries_sample) {
    int numSamples = (int)(num_entries * random_entries_sample_factor);
    sequence<int> entriesSampleIdxs = tabulate(numSamples, [&] (size_t i) {
      thread_local std::random_device r;
      thread_local std::default_random_engine e1(r());
      thread_local std::uniform_int_distribution<int> rng(0, num_entries - 1);
      return rng(e1) % num_entries;
    });

    B = map(A, [&] (feature f) {
      row v = tabulate(numSamples, [&] (size_t i){
        return f.vals[entriesSampleIdxs[i]];
      });
      value maxv = reduce(v, maxm<value>());
      feature fp = feature(f.discrete, maxv + 1, v);
      return fp;
    }, 1);
  }
  std::set<int, std::greater<int>> usedFeatures;
  return build_tree(B, verbose, usedFeatures);
}

auto buildForest(features const &Train, bool verbose) {
  return tabulate(num_trees, [&] (size_t i) {
    return buildRandomTree(Train, verbose);
  }, 1);
}

row classifyRandomForest(features const &Train, rows const &Test, bool verbose) {
  sequence<value> testLabels = map(Test, [&] (row const &r) {
    return r[0];
  });
  value maxLabel = std::max(Train[0].num, (int)reduce(testLabels, maxm<value>()));

  sequence<tree*> forest = buildForest(Train, verbose);

  return map(Test, [&] (row const &r) {
    sequence<value> preds = map(forest, [&] (tree* T) -> value { return classify_row(T, r); }, 1);
    value maj = majority(preds, maxLabel);
    return maj;
  }); 
}

row classify(features const &Train, rows const &Test, bool verbose) {
  auto startTime = time(0);
  row result;
  if(use_random_forest) {
    result = classifyRandomForest(Train, Test, verbose);
  }
  else {
    features A = Train;
    std::set<int, std::greater<int>> usedFeatures;
    tree* T = build_tree(A, verbose, usedFeatures);
    if (true) cout << "Tree size = " << T->size << endl;
    int num_features = Test[0].size();
    result = map(Test, [&] (row const& r) -> value {return classify_row(T, r);});
  }
  auto endTime = time(0);
  cout << "Time: " << endTime - startTime << endl;
  return result;
}
