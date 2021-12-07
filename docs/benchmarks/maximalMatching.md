---
title: Maximal Matching
---

# Maximal Matching (MM):

Given a undirected graph return a maximal matching for the graph.
The input graph can be in any format (as long as it does not
encode the MIS somehow).  Also the code cannot reorder the graph for
locality.  The output needs to be a sequence of integers corresponding
to the positions of the edges in the input (zero based) in the MM.
The output sequence need not be sorted.

###  Input Distributions

The default input graphs are as follows:

- Random local graph with approximately n vertices and 20n
edges.   It should be generated using:  
`randLocalGraph -d 3 -m <10*n> <n> <filename>`  
`n` = 20 million for large instances and 2 million for small.

- An RMAT graph with parameters (.2,.125,.125,.55) and about 12n edges.
It should be generated using:  
`rMatGraph -a .55 -b .125 -m <12*n> <n> <filename>`  
`n` = 16 million for large instances and 2.25 million for small.

- A 3 dimensional grid generated using:  
`gridGraph -d 3 <n> <filename>`  
`n` = 64 million for large instances and 8 million for small.

### Input and Output File Formats 

The input is a graph in the in the edge graph format.  The 
output needs to be in the sequence format.  
