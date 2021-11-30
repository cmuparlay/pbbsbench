---
title: The PBBS Benchmark Suite (V2)
---

New version (2020) of pbbs benchmarks

This repository uses a submodule (parlaylib).    To pull the submodule, in pbbsbench:

\> git submodule init

\> git submodule update

To run a test, for example, try:

\> cd comparisonSort/sampleSort

\> make

\> ./testInputs -r 3

Where -r 3 means run three times

The other options
   [-p n] : run on n processors
   [-x] : do not check the result

[Comparison Sort]({% link foo.md %})
