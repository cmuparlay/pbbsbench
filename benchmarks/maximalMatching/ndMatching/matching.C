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

#define NOTMAIN 1
#include <iostream>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/graph.h"
#include "matching.h"

parlay::sequence<edgeId> maximalMatching(edges const &E) {
  size_t n = std::max(E.numCols,E.numRows);
  size_t m = E.nonZeros;
  enum state : char {Init, Locked, Taken};
  auto matched = parlay::tabulate<std::atomic<state>>(n, [] (int i) {return Init;});
  parlay::sequence<bool> result(m, false);
  parlay::parallel_for(0, m, [&] (int i) {
        long u = E[i].u;
        long v = E[i].v;
        if (u == v) return;
        while (true) {
          state old = Init;
          if (matched[u].load() == Taken || matched[v].load() == Taken) {
            result[i] = false;
            return;
          } else if (matched[u].compare_exchange_strong(old, Locked)) {
            if (matched[v].compare_exchange_strong(old, Taken)) {
              matched[u] = Taken;
              result[i] = true;
              return;
            } else matched[u] = Init;
          }
        }});
  auto matchingIdx = parlay::pack(parlay::delayed::tabulate(m, [&] (unsigned int i) {return i;}), result);
  return matchingIdx;
}  
