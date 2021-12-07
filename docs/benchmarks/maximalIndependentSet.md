---
title: Maximal Independent Set
---

# Maximal Independent Set (MIS)

Given a undirected graph return a maximal independent set for the
graph.  The input graph can be in any format (as long as it does not
encode the MIS somehow).  Also the code cannot reorder the graph for
locality.  The output needs to be a sequence of the vertices appearing
in the MIS.  It need not be in sorted order.

<h3>Default Input Distributions</h3>

The default input graphs are as follows:

- Random local graph with approximately n vertices and 20n
edges.   It should be generated using:  
`randLocalGraph -j -d 3 -m <10n> <n> <filename>`  
`n` = 20 million for large instances and 2 million for small.

- An RMAT graph with parameters (.2,.125,.125,.55) and about 12n edges.
It should be generated using:  
`rMatGraph -j -a .55 -b .125 -m <12n> <n> <filename>`  
`n` = 16 million for large instances and 2.25 million for small.

- A 3 dimensional grid generated using:  
`gridGraph -j -d 3 <n> <filename>`  
`n` = 64 million for large instances and 8 million for small.

### Input and Output File Formats

The input is a graph in the in the adjacency graph file format.  The
output is a sequence of integers corresponding to the locations of all
vertices in the MIS (0 based).    The output needs to be in the sequence format. 

