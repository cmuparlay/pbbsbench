// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2010 Guy Blelloch and the PBBS team
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
#include <cstring>
#include "parlay/parallel.h"
#include "parlay/internal/quicksort.h"
#include "common/sequenceIO.h"
#include "common/parseCommandLine.h"
using namespace std;
using namespace benchIO;

template <typename T, typename LESS, typename Key>
void check_sort(sequence<sequence<char>> In,
		sequence<sequence<char>> Out,
		LESS less, Key f) {
  sequence<T> in_vals = parseElements<T>(In.cut(1, In.size()));
  sequence<T> out_vals = parseElements<T>(Out.cut(1, In.size()));
  size_t n = in_vals.size();
  auto sorted_in = parlay::stable_sort(in_vals, less);
  parlay::internal::quicksort(make_slice(in_vals), less);

  atomic<size_t> error = n;
  parlay::parallel_for (0, n, [&] (size_t i) {
    if (f(in_vals[i]) != f(out_vals[i])) 
	parlay::write_min(&error,i,std::less<size_t>());
  });
  
  if (error < n) {
    auto expected = parlay::to_chars(f(in_vals[error]));
    auto got = parlay::to_chars(f(out_vals[error]));
    cout << "comparison sort: check failed at location i=" << error
	 << " expected " << expected << " got " << got << endl;
    abort();
  }
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"<infile> <outfile>");
  pair<char*,char*> fnames = P.IOFileNames();
  char* infile = fnames.first;
  char* outfile = fnames.second;

  auto In = get_tokens(infile);
  elementType in_type = elementTypeFromHeader(In[0]);
  size_t in_n = In.size() - 1;

  auto Out = get_tokens(outfile);
  elementType out_type = elementTypeFromHeader(Out[0]);
  size_t out_n = In.size() - 1;

  if (in_type != out_type) {
    cout << "sortCheck: types don't match" << endl;
    return(1);
  }
  if (in_n != out_n) {
    cout << "sortCheck: lengths dont' match" << endl;
    return(1);
  }

  if (in_type == doubleT) {
    check_sort<double>(In, Out, std::less<double>(), [&] (double x) {return x;});
  } else if (in_type == doublePairT) {
    using dpair = pair<double,double>;
    auto less = [] (dpair a, dpair b) {return a.first < b.first;};
    check_sort<dpair>(In, Out, less, [&] (dpair x) {return x.first;});
  } else if (in_type == stringT) {
    using str = sequence<char>;
    auto strless = [&] (str const &a, str const &b) {
      auto sa = a.begin();
      auto sb = b.begin();
      auto ea = sa + min(a.size(),b.size());
      while (sa < ea && *sa == *sb) {sa++; sb++;}
      return sa == ea ? (a.size() < b.size()) : *sa < *sb;
    };
    check_sort<str>(In, Out, strless, [&] (str x) {return x;});
  } else if (in_type == intType) {
    check_sort<int>(In, Out, std::less<int>(), [&] (int x) {return x;});
  } else {
    cout << "sortCheck: input files not of accepted type" << endl;
    return(1);
  }
}
