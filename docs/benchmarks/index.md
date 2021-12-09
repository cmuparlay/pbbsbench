---
title: PBBS Listing of Benchmarks 
---

#  PBBS Listing of Benchmarks

Here we include a listing for each benchmark.   Follow the link for
the benchmark to get more details, including the input/output
specification and the default distributions.

We break the benchmarks into four broad categories, and other.

### Basic Building Blocks

- [comparisonSort](comparisonSort.html) (SORT)  
Returns the sorted input based on a comparison-based sort.

- [histogram](histogram.html) (HIST)  
Returns the histogram for a sequence of integers.

- [integerSort](integerSort.html) (ISORT)  
Sorts a sequence of integers, possibly with tag-along values. 

- [removeDuplicates](removeDuplicates) (DDUP)  
Returns the input sequence with duplicates removed.

### Graph Algorithms

- [breadthFirstSearch](breadthFirstSearch.html) (BFS)  
Returns a breadth-first-search tree from a given vertex in a graph.

- [maximalIndependentSet](maximalIndependentSet.html) (MIS)  
Returns a maximal independent set for an undirected graph. 

- [maximalMatching](maximalMatching.html) (MM)  
Returns a maximal matching for an undirected graph. 

- [minSpanningForest](minSpanningForest.html) (MSF)  
Returns a minimum spanning forest of a weighed undirected graph. 

- [spanningForest](spanningForest.html) (SF)  
Returns a spanning forest of an undirected graph. 

### Text Processing

- [BWDecode](BWDecode.html) (BWD)  
Decodes a string encoded with the Burrows-Wheeler transform.

- [invertedIndex](invertedIndex.html) (IIDX)  
Returns an inverted index given  a string of documents.

- [longestRepeatedSubstring](longestRepeatedSubstring.html) (LRS)  
Returns the longest repeated substring in a string.

- [suffixArray](suffixArray.html) (SA)  
Returns the suffix array for a string. 

- [wordCounts](wordCounts.html) (WC)  
Counts the number of occurrences of each word in a string. 

### Computational Geometry/Graphics

- [convexHull](convexHull.html) (CH)  
Returns the convex hull for a set of points in 2 dimensions.

- [delaunayRefine](delaunayRefine.html) (DR)  
Adds points to a Delaunay triangulation (DT) in 2d, so resulting DT has no 
small angles. 

- [delaunayTriangulation](delaunayTriangulation.html) (DT)  
Returns the Delaunay triangulation of points in 2d. 

- [nearestNeighbors](nearestNeighbors.html) (KNN)  
Returns the k nearest neighbors for points in 2d and 3d. 

- [rayCast](rayCast.html) (RAY)  
For a set of rays and a set of triangles returns the first triangle 
each ray intersects. 

- [rangeQuery2d](rangeQuery2d) (RQ)  
For a set of points and a set of rectangles, returns a count of the number of 
points in each rectangle. 

### Other

- [classify](classify.html) (CLAS)  
Predicts labels for feature vectors based on a training vectors with attached labels.

- [nBody](nBody.html) (NBODY)  
Returns gravitational forces among points in 3d. 
