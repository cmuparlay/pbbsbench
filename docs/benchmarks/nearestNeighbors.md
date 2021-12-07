---
title: K-Nearest Neighbors
---

# K-Nearest Neighbors (KNN)

Given n points in 2 and 3 dimensions find the k nearest neighbors for
each point based on Euclidean distance.  The input is a sequence of
points, each of which is a pair (2d) or tripple (3d) of
double-precision floating-point numbers.  The value k is specified as
a parameter to the program.  The output is a sequence n tuples of
length k each.   Each tuple identifies the indices of its k nearest
neighbors from the input sequence (zero based).

### Default Input Distributions

The default distributions are the following:

- k = 1. Points in 2 dimensions selected uniformly at random in 
  a square.   Should be generated with__
`randPoints -d 2 <n> <filename>`

- k = 1. Points chosen at random in 2 dimensions from the Kuzmin distribution.   Should be 
generated with:  
`randPoints -k -d 2 <n> <filename>`.

- k = 1. Points in 3 dimensions selected uniformly at random in a cube. 
Should be generated with:  
`randPoints -d 3 <n> <filename>`. 

- k = 1. Points chosen at random in 3 dimensions from the surface of a 
sphere.   Should be generated with:  
`randPoints -S -d 3 <n> <filename>`.

- k = 10. Points in 3 dimensions selected uniformly at random in 
  a cube.   Should be generated with  
`randPoints -d 3 <n> <filename>`

- k = 10. Points chosen at random in 3 dimensions from the Plummer distribution.   Should be 
generated with:  
`randPoints -k -d 3 <n> <filename>`.

The large size is n = 10 million, and the small size is n = 1 million.

# Input and Output File Formats

The input needs to be in the [points file format](../fileFormats/geometry.html#points).
The output needs to be in the [sequence file
format](../fileFormats/sequence.html) with a total of (n x k)
entries and integer type.  Each point corresponds to k adjacent
entries, where points are ordered as in the input.
