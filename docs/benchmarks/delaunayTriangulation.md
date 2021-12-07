---
title: Delaunay Triangulation
---

# Delaunay Triangulation (DT)

Given a set of points in 2 dimensions generate the Delaunay
triangulation.  The input should be a sequence of points, each a pair
of double precision floating-point numbers.  By default one can assume
the points are in general position (no three on a line or four on a
circle) and that double precision arithmetic for in-circle tests
suffice.  Our input sets satisfy these conditions.  The purpose of
this benchmark is not to compare state-of-the-art exact geometric
predicates.  If the implementation is robust in general, then this can
be noted.

The triangulation can add points.  Such points can be used, for
example, to insert points at the corners of a bounding triangle around the
original points.   The triangulation needs to be a proper Delaunay
triangulation of all points.

The output must be a sequence of points and a sequence of triangles.
The points should start with the original points and can have the
extra points at the end.  The triangles should be represented as
tripple of integer indices indicating the position of its three corner
points in the input sequence (zero based).  Triangles can be in any
order.

### Default Input Distributions

The distributions are:

- Points chosen uniformly at random within a unit circle.   Should be 
generated with:  
`randPoints -s -d 2 <n> <filename>`.

- Points chosen at random from the Kuzmin distribution.   Should be 
generated with:  
`randPoints -k -d 2 <n> <filename>`.   

The large size is n = 10 million, and the small size is n = 1 million.

### Input and Output File Formats

The input needs to be in the [2dpoints file format](../fileFormats/geometry.html#points).
The output needs to be in [triangle file format](../fileFormats/geometry.html#triangles).
