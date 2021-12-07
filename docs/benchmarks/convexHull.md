---
title: Convex Hull
---

# Convex Hull (CH)

Given a set of points in 2 dimensions return the set of points on the
convex hull.  The output is the index of the points on the hull in
clockwise sorted order around the hull starting with the point with
minimum x-coordinate.  Each point must be a pair of double-precision
floating-point numbers.  The points are generated so it is unlikely
that a double-precision line-side test will be ambiguous due to
roundoff (i.e., exact arithmetic is not required).

### Default Input Distributions

The default distributions are the following:

- Points chosen uniformly at random within a unit circle.   Should be 
generated with:  
`randPoints -s -d 2 <n> <filename>`.

- Points chosen at random from the Kuzmin distribution.   Should be 
generated with:  
`randPoints -k -d 2 <n> <filename>`.

- Points chosen uniformly at random on the perimeter of a unit circle.   Should be 
generated with:  
`randPoints -S -d 2 <n> <filename>`.

The large size is n = 100 million, and the small size is n = 10 million.

### Input and Output File Formats

The input needs to be in the [2dpoints file format](../fileFormats/geometry.html#points).
The output needs to be in the [sequence file format](../fileFormats/sequence.html)
