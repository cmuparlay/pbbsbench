---
title: Integer Sort
---

# Integer Sort (ISORT)

Sort fixed-length unsigned integer keys into ascending order with the ability to carry
along fixed-length auxiliary data.  The sort must be stable unless 
the benchmark is indicated with -U (ISORT-U).    The number of bits to
sort can be specified to the program.

### Default Input Distributions

The default distributions are as follows:

- n unsigned integers generated uniformly at random in the range [0:n).
Should be generated with:  
`randomSeq -t int <n> <filename>`

- n unsigned integers generated from the exponential distribution in the 
range [0:n).
Should be generated with:  
`exptSeq -t int <n> <filename>`

- n unsigned integers generated uniformly at random in the range 
  [0:n) each tagged with data also in the range [0:n). 
  Should be generated with:  
  `randomSeq -t int <n> <tmpfile>`  
  `addDataSeq -t int <tmpfile> <filename>`

- n unsigned integers generated uniformly at random in the range 
  [0:256) each tagged with data in the range [0:n). 
  Should be generated with:  
  `randomSeq -t int -r 256 <n> <tmpfile>`  
  `addDataSeq -t int <tmpfile> <filename>`

The large size is n = 100 million, and the small size is n = 10
million.

### Input and Output File Formats

The input and output data need to be in the [sequence file format](../fileFormats/sequence.html),
both with the same element type. The element type is either a pair of
integers or a single integer.

The output file must be in sorted order with respect to integer
ordering (first integer if pairs).
