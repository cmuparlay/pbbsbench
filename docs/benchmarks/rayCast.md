---
title: Ray-Triangle Intersection
---

# Ray-Triangle Intersection (RAY)

Given a set of triangles contained inside a 3d bounding box and a set
of rays that penetrate the box, calculate for each ray the first
triangle it intersects, if any.

The input consists of a sequence of triangles and a sequence of rays.
The triangles are defined by a sequence of points (each a tripple of
double-precision floating-point numbers) and a sequence of integer
tripples.  Each integer tripple represents the index of its three
corners in the point sequence (zero based).  All points are within a
bounding box.  The rays are each defined as a point in 3d along with a
vector in 3d.  The ray starts at the point and points in the direction
of the vector.  All 6 values use double-precision floating-point.  All
rays start at or to the left of the bounding box (by x-coordinate),
and point right.

The output is a sequence corresponding to the rays (in the same order)
along with the index of the first triangle each ray intersects (using
the ordering of the triangles, zero based).   If the ray does not
intersect a triangle then the index should be -1.

### Default Input Distributions

This benchmark is run on three meshes, which are distributed as part
of the benchmark.  In all cases a set of random rays coming from the
left side of the bounding box are generated using

```
addRays <triangelFile> <rayFile>
```
This generates a number of rays equal to the number of triangles.
The input files are:

  - Happy.  This input is the famous "happy Buddha" mesh from the
  [Stanford 3D Scanning
  Repository](http://graphics.stanford.edu/data/3Dscanrep/).  
  Consists of 1087716 triangles.

  - Angel. This input is from Georgia Tech.  
  Consists of 474048 triangles.

  - Dragon.   This input is also from the
  [Stanford 3D Scanning
  Repository](http://graphics.stanford.edu/data/3Dscanrep/).  
  Consists of 871414 triangles.


### Input and Output File Formats

The input consist of two files.   The first is in [triangle file
  format](../fileFormats/geometry.html#triangles), and the second in
  the [3dpoints file format](../fileFormats/geometry.html#points).
  For n points, the second file contains 2n entries where adjacent
  pairs correspond to the point-vector pair for each ray.
  The output needs to be in the [sequence file
  format](../fileFormats/sequence.html)
  with integer types.
