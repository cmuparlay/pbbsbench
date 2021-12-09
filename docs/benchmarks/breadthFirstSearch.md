---
title: Breadth First Search
---

# Breadth First Search (BFS)

Run a breadth first search on a directed graph starting at a specified
source vertex and return a bfs-tree.  Duplicate input edges are
allowed.  The input is in adjacency array format with an ordering
among the outgoing edges.  If the graph is not connected then only the
tree over the vertices reachable from the source should be returned.

The input graph can be represented as desired (probably some form of
adjacenty list or adjacency array).   The graph
can be directed.   The search will start at a specified start vertex
and search everything reachable from the start.

For a graph with n verticels the output is an integer sequence of
lenght n representing a valid BFS tree from the start.  In particular
each location will be the index of its parent in a BFS tree resulting
from the search.  The root will point to itself, and any unvisted
vertices will have -1 as their entry.

### Default Input Distributions

The default input graphs are as follows:

- Random local graph with approximately n vertices and 20* n
edges.   It should be generated using:  
`randLocalGraph -j -d 3 -m <10n> <n> <filename>`  
`n` = 20 million for large instances and 2 million for small.

- An RMAT graph with parameters (.2,.125,.125,.55) and about 12*n edges.
It should be generated using:  
`rMatGraph -j -a .55 -b .125 -m <12n> <n> <filename>`  
`n` = 16 million for large instances and 2.25 million for small.

- A 3 dimensional grid generated using:  
`gridGraph -j -d 3 <n> <filename>`  
`n` = 64 million for large instances and 8 million for small.

### Input and Output File Formats

The input is a graph in the in the [adjacency graph format](../fileFormats/graph.html)
The output will be in be an [integer sequence format](../fileFormats/sequence.html)
