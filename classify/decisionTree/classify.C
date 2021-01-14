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
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/internal/get_time.h"
#include "classify.h"

using namespace parlay;
double infinity = std::numeric_limits<double>::infinity();

struct tree {
  bool is_leaf;
  int feature_index;
  int feature_cut;
  value best;
  sequence<tree*> children;
  tree(int i, int c, sequence<tree*> children)
    : is_leaf(false), feature_index(i), feature_cut(c), children(children) {}
  tree(value best) : is_leaf(true), best(best) {}
};

tree* Leaf(int best) {
  return new tree(best);
}

tree* Node(int i, int cut, sequence<tree*> children) {
  return new tree(i, cut, children);
}

// To be put into parlay

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
  auto x = histogram(a,m);
  return max_element(x) - x.begin();
}

// End to be put into parlay

// entropy * total given a histogram a (total = sum(a))
template <typename Seq>
double entropy(Seq a, int total) {
  return reduce(map(a, [=] (int l) {
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
    double e = entropy(low_counts, m) + entropy(high_counts, n - m)
      + log2(float(i+1)) + log2(float(b.num-i-1));
    if (e < cur_e) {
      cur_e = e;
      cur_i = i+1;
    }
  }
  return std::pair(cur_e, cur_i);
}

// information content of s (i.e. entropy * size)
double info(values s, int num_vals) {
  size_t n = s.size();
  if (n == 0) return 0.0;
  auto x = histogram(s, num_vals);
  return entropy(x, n);
}

// info of a conditioned on b
double cond_info_discrete_o(feature const &a, feature const &b) {
  auto groups = group_by_index(delayed_zip(b.vals, a.vals), b.num);
  return reduce(map(groups, [&] (auto g) {return log2(float(1 + g.size())) + info(g, a.num);}));
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

auto build_tree(features &A) {
  int num_features = A.size();
  int num_entries = A[0].vals.size();
  if (num_entries < 64 || all_equal(A[0].vals))
    return Leaf(majority(A[0].vals, A[0].num));

  double label_info = info(A[0].vals,A[0].num);
  auto costs = tabulate(num_features - 1, [&] (int i) {
      if (A[i+1].discrete) {
	return std::tuple(cond_info_discrete(A[0], A[i+1]), i+1, -1);
      } else {
	auto [info, cut] = cond_info_continuous(A[0], A[i+1]);
	return std::tuple(info, i+1, cut);
      }},1);

  auto min1 = [&] (auto a, auto b) {return (std::get<0>(a) < std::get<0>(b)) ? a : b;};
  auto min_m = make_monoid(min1, std::tuple(infinity, 0, 0));
  auto [best_info, best_i, cut] = reduce(costs, min_m);
  double threshold = log2(float(num_features));

  //cout << num_entries << ", " << best_i << ", " << cut << ", " << label_info << ", " 
  //     << best_info << endl;

  if (label_info - best_info < threshold)
    return Leaf(majority(A[0].vals, A[0].num));
  else {
    int m;
    values split_on;
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

    auto r = map(B, [&] (features &a) {return build_tree(a);}, 1);
    return Node(best_i, cut, r);
  }
}

void classify(features &A) {
  build_tree(A);
}
