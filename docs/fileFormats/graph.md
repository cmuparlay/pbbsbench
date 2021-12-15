---
title: Graph File Formats
---

# Graph File Formats

The graph file formats include an **adjacency graph** format, an
**edge graph** format, and a **weighted edge graph** format.  All
formats are ascii and entries are delimited by any consecutive
sequence of delimiter characters: **tab**, **space**, **line feed**
(ascii 0x0A), and **carriage return** (ascii 0x0D).  Files can start
and end with delimiters, which are ignored.  Throughout the
description n refers to the number of vertices and m to the number of
edges in a graph.  

### Adjacency Graph

The adjacency graph format starts with a sequence of offsets one for
each vertex, followed by a sequence of directed edges ordered by their
source vertex.  The offset for a vertex i refers to the location of
the start of a contiguous block of out edges for vertex i in the
sequence of edges.  The block continues until the offset of the next
vertex, or the end if i is the last vertex.  All vertices and offsets
are 0 based and represented in decimal.  The specific format is as
follows:

```
AdjacencyGraph
<n>
<m>
o0
o1
...
o(n-1)
e0
e1
...
e(m-1)
```

### Edge Graph

The edge graph format consists of a sequence of edges/arcs each being
a pair of integers. 
The format can either be interpreted as a directed graph or an undirected graphs
depending on the application.
Vertices are assumed to start at 0.
The specific format is as follows:

```
EdgeArray
s0 t0 
s1 t1 
... 
s(n-1) t(n-1) 
```

where `si` and `ti`
refer to the two endpoints of the undirected edge i, or the source and
target of a directed edge (arc) i.

### Weighted Edge Graph

The weighted edge graph format is the same as the edge graph format but
includes double-precision floating-point weights.
The specific format is as follows:

```
WeightedEdgeArray
s0 t0 w0
s1 t1 w1
... 
s(n-1) t(n-1) w(n-1)
```

where `wi` is the weight of edge i.  The weight can either
be in decimal or exponential notation.
