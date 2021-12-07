---
title: Delaunay Refinement
---

# Delaunay Refinement (DR)

Given a Delaunay triangulation in 2 dimensions, return a Delaunay triangulation
that includes all original points plus additional points such that no triangle
has an angle smaller than a threshold **Theta**.
The input is a sequence of points, each a pair of double-precision
floating-point numbers, and a sequence of triangles.   The triangles
are a tripple of indices, indicating the position of their corners in
the point sequence (zero-based).   The input must be a valid Delaunay
triangulation.

The output is another valid Delaunay triangulation in the same
format.  The points must start with the input points in the same
order, and can contain additional points after them.   Furthermore
the output triangles cannot have any angles smaller than
**Theta**.

## Default Input Distributions

For this benchmark the input triangulation is generated using a
Delaunay triangulation on various point distributions.  These
triangulations are generated using the PBBSbench code for Delaunay
triangulations, which also adds in a small number of boundary points
beyond the original input points.  The point sets are non degenerate.

- Points chosen uniformly at random within a unit circle.   Should be 
generated with  
`randPoints -s -d 2 <n> <tmpfile>`  
`delaunay -o <filename> <tmpfile>`

- Points chosen at random from the Kuzmin distribution.   Should be 
generated with  
`randPoints -k -d 2 <n> <tmpfile>`  
`delaunay -o <filename> <tmpfile>`

The large size is n = 10 million, and the small size is n = 1 million.

### Input and Output File Formats

The input and output need to be in [triangle file format](../fileFormats/geometry.html#triangles).
