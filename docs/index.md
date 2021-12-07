---
title: The PBBS Benchmark Suite (V2)
---

#  The PBBS Benchmark Suite (V2)

This documents the Problem Based Benchmark Suite (PBBS), a collection
of over 20 benchmarks defined in terms of their IO characteristics.
They are designed to make it possible to compare different algorithms,
or implementations in different programming languages.

For each benchmark the suite provides:

- the speficication of the input and expected output for the problem
- the specification of a default set of input instances 
- code for generating inputs (written to a file)
- code for checking correctness of output (read from a file)
- code for timing the benchmark across the instances
- a default parallel implementation
- a default sequential implementation (for most benchmarks)
- a variety of other implementations (for some benchmarks)

The benchmarks are designed to be agnostic to the programming
language.  However, the framework is mostly in C++ and some of the
tools are easier to use with C++.  For example there is a timing
driver for C++ that can be linked with any implementation of a
benchmark.

## Getting Started

The benchmark suite can be downloaded from github using:

```
> git clone https://github.com/cmuparlay/pbbsbench.git
```

It uses two submodules, which can be initialized with:

```
> git submodule init
> git submodule update
```

### Requirements

The benchmarks have been tested on Ubuntu and MacOS.

The software requirements are:

- C++-17 compiler (tested with gcc and clang)

The system requirements are

- for small data (12B of RAM and 10GB of disk)
- for large data (64GB of RAM and 90GB of disk)

The following are not required, but will give better performance

- jemalloc  (only gives slight performance improvement)
- 20+ cores (the more cores the faster)
- numactl installed (if this is not installed you need to run "./runall -nonuma")

### Running the benchmarks

The command `./runall` will compile and run all the benchmarks
reported but will take a couple hours.  For a faster run, try:

```
  ./runall -par -small
```
  
This only compiles and runs the parallel benchmarks, which run much faster, and on
significantly smaller input data (an order of magnitude smaller for some benchmarks).
This runs in 12 minutes on a 20 core (40 hyperthread) machine.

You can also test individual benchmarks.   For example, you can test the
parallel comparison sort using:

```
  ./runall -only comparisonSort/sampleSort
 ```
  
This will run the parallel sampleSort on the default (full sized) inputs.
The call

```
  ./runall -small -only comparisonSort/sampleSort
```
  
will run it on the smaller inputs.  More details on arguments are
given below.

## More Details


### Options

The ./runall has the following options which can be extracted by using
`./runall -h`.

```
  -scale    : this runs it on a range of different thread counts up the the number of threads on the machine
  -small    : runs tests on smaller inputs (calls ./testInput_small instead of ./testInput).
  -par      : only run benchmarks that are parallel (saves time)
  -only <name>   : only run a particular benchmark
  -notime   : only compile the benchmarks
  -nonuma   : don't use numactl
  -nocheck  : don't check correctness of results (saves time)
```
  
For the `-only` option use the path to the implementation, e.g.

```
  ./runall -only comparisonSort/sampleSort
```

### Directory Organization

The benchmarks themselves are organized hierarchically.  At the top
level are the following benchmarks:

- [breadthFirstSearch](benchmarks/breadthFirstSearch.html) (BFS) --
  Returns a breadth-first-search tree from a given vertex in a graph.
- [BWDecode](benchmarks/BWDecode.html) (BWD) -- Decodes a string
  encoded with the Burrows-Wheeler transform.
- classify (CLAS) -- Predicts labels for a test set of feature vectors
  based on a training set with attached labels.
- [comparisonSort](benchmarks/comparisonSort.html) (SORT) --
  Returns the sorted input based on a comparison-based sort.
- [convexHull](benchmarks/convexHull.html) (CH) -- Returns the
  convex hull for a set of points in 2 dimensions.
- [delaunayRefine](benchmarks/delaunayRefine.html) (DR) -- Given a
  Delaunay triangulation (DT) in 2d, adds points so resulting DT has no
  small angles.
- [delaunayTriangulation](benchmarks/delaunayTriangulation.html) (DT)
  -- Returns the Delaunay triangulation of points in 2d.
- [histogram](benchmarks/histogram.html) (HIST) -- Returns the
  histogram for a sequence of integers.
- [index](benchmarks/invertedIndex.html) (II) -- Returns an inverted
  index given  a string of documents.
- [integerSort](benchmarks/integerSort.html) (ISORT) -- Sorts a
  sequence of integers, possibly with tag-along values.
- longestRepeatedSubstring (LRS) -- Returns the longest repeated
  substring in a string.
- [maximalIndependentSet](benchmarks/maximalIndependentSet.html) (MIS)
  -- Returns a maximal independent set for an undirected graph.
- [maximalMatching](benchmarks/maximalMatching.html) (MM) -- Returns a
  maximal matching for an undirected graph.
- [minSpanningForest](benchmarks/minSpanningForest.html) (MSF) --
Returns a minimum spanning forest of a weighed undirected graph.
- [nBody](benchmarks/nBody.html) (NBODY) -- 
Returns gravitational forces among points in 3d. 
- [nearestNeighbors](benchmarks/nearestNeighbors.html) (KNN) --
Returns the k nearest neighbors for points in 2d and 3d.
- [rayCast](benchmarks/rayCast.html) (RAY) --
For a set of rays and a set of triangles returns the first triangle
each ray intersects.
- rangeQuery2d (RQ) --
For a set of points and a set of rectangles, returns a count of the number of
points in each rectangle.
- [removeDuplicates](benchmarks/removeDuplicates) (DDUP) --- Returns the input sequence with
  duplicates removed.
- [spanningForest](benchmarks/spanningForest.html) (SF) -- Returns a
  spanning forest of an undirected graph.
- [suffixArray](benchmarks/suffixArray.html) (SA) -- Returns the suffix array for a string.
- [wordCounts](benchmarks/wordCounts.html) (WC) -- Counts the number of occurrences of each word in a string.

Within each benchmarks are the implementations.   Each benchmarks also
have some directories shared across implementations.  In particular
each one has a directory called `bench` which contains the driver
code, testing code, and specification of default inputs.  There is
also a directory for generating the data.

Within each implementation directory, you can run `make` to make the
executable, and then run `./testInputs` to run the benchmarks.  These
are run automatically by the ./runall script.  On a machine with
multiple chips, using `numactl -i all ./testInputs` will give better
results.  `./testInputs_small` will use the smaller inputs.

The `testInputs` script has several options including:

```
  -x : do not check the output
  -r <count>  : number of rounds to use
  -p <count>  : number of threads to use
  ```
  
The actual inputs are specified in the script and can be changed if desired.



