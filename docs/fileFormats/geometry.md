---
title: Geometry File Formats
---

# Geometry File Formats

The geometry file formats include **points** in 2 and 3
dimensions and **triangles**.  All formats are ascii and
entries are delimited by any consecutive sequence of delimiter
characters: **tab**, **space**, **line
feed** (ascii 0x0A), and **carriage return**
(ascii 0x0D).  Files can start and end with delimiters, which are
ignored.

### Points

The points file format consists of a sequence of pairs of coordinates
(2 dimensions) or triples of coordinates (3 dimensions).  Each
coordinate is a real number capable of fitting withing double
precision accuracy.  The number can be either in decimal or
exponential notation.
The specific format for 2d points is:

```
pbbs_sequencePoint2d 
<x0> <y0>
<x1> <y1>
...
<x_(n-1)> <y_(n-1)>
```

and for 3d points:

```
pbbs_sequencePoint3d
<x0> <y0> <z0>
<x1> <y1> <z1>
...
<x_(n-1)> <y_(n-1)> <z_(n-1)>
```

### Ray

For rays we use the same format as points, but pair up adjacent points.
The first of each pair is the source of the ray and the second the direction.
Therefore a file with n rays is equivalent to one with 2n points.

### Triangles

The triangles format starts with a sequence of points and is
followed by a sequence of triangles.  Each point is either in 2d or 3d
and stored as in the **points** format.  Each triangle
is a triple consisting of the indices of the points at its three
corners.  The indices are zero based and refer to the position of the
point in the file.  The number of points and triangles is
included in the file.  In particular for 2d points the format is
as follows:

```
pbbs_triangles
<n>
<m>
<x0> <y0>
<x0> <y0>
<x1> <y1>
...
<x_(n-1)> <y_(n-1)>
<a0> <b0> <c0>
<a1> <b1> <c1>
...
<a_(m-1)> <b_(m-1)> <c_(m-1)>
```

For 3d points each point has an an additional `<zi>`.
