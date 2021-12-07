---
title: Minimum Spanning Forest
---

# Minimum Spanning Forest (MSF)

Given a weighted undirected graph return the minimum spanning tree
(MST), or minimum spanning forest (MSF) if the graph is not connected.
The input graph can be in any format (as long as it does not encode
the MIS somehow).  Also the code cannot reorder the graph for
locality.  The weights must be single-precision floating-point
numbers.  The output needs to be a sequence of integers corresponding
to the positions of the edges in the input (zero based) that are in
the ST.  The output sequence need not be sorted.  In the case that
there are multiple possible MSFs (due to equal weights), any MSF is
valid.

### Default Input Distributions

The default input graphs are as follows:

- Random local graph with approximately n vertices and 20* n
edges.   It should be generated using:  
`randLocalGraph -d 3 -m <10*n> <n> <tmpfile>`  
`addWeights <tmpfile> <filename>`  
`n` = 20 million for large instances and 2 million for small.

- An RMAT graph with parameters (.2,.125,.125,.55) and about 12*n edges.
It should be generated using:  
`rMatGraph -a .55 -b .125 -m <12*n> <n> <tmpfile>`  
`addWeights <tmpfile> <filename>`  
`n` = 16 million for large instances and 2.25 million for small.

- A 3 dimensional grid generated using:  
`gridGraph -d 3 <n> <tmpfile>`  
`addWeights <tmpfile> <filename>`  
`n` = 64 million for large instances and 8 million for small.

### Input and Output File Formats

The input is a graph in the in the weighted edge graph file format. 
The output needs to be in the sequence file format.

