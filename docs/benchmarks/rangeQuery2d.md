---
title: 2d Range Query
---

# 2D Range Query (RAY)

Given a set of points and a set of rectangles in 2 dimensions, return
for each rectangle  a count of how many points are contained in the
rectangle.

The input consists of a sequence of n points, each a pair of
double-precision floating point numbers, and a sequence of m rectangles,
each consisting of two points (i.e. two double-precision fp numbers
each).  The two points for a rectangle represent the lower-left corner
and the upper-right corner.

The output is a sequence of m integers corresponding to the count for
each rectangle.   Any point on the boundary of a rectangle counts as
inside the rectangle.

### Default Input Distributions

The default distributions use a total of n points, where the first 2
(n / 3) points are used for the query rectangles, where each rectangle
uses two adjacent points.  Here we assume (n / 3) rounds down.  The
rest of the points (n - 2 (n/3) of them) are used as points to query
over.

The default distributions are the following:

- Points chosen uniformly at random within a unit cube.
Should be 
generated with:  
`randPoints -d 2 <n> <filename>`.

- Points chosen at random from the Kuzmin distribution.   Should be 
generated with:  
`randPoints -k -d 2 <n> <filename>`.

The large size is n = 10 million, and the small size is n = 1 million.


### Input and Output File Formats


The input needs to be in the [2dpoints file format](../fileFormats/geometry.html#points).
The output needs to be in the
[sequence file format](../fileFormats/sequence.html) with integer type.
