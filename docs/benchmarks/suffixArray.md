---
title: Suffix Arrays
---

# Suffix Arrays (SA)

Given a string generate its [suffix
array](http://en.wikipedia.org/wiki/Suffix_array), i.e. the
sorted sequence of all suffixes of the input.

The input is a string of length n containing no null characters, and
the output is the suffix array as a sequence of length n.    The
indices in the sequence are zero-based (i.e. the location of 0 in the
array gives the rank of the whole string among all its suffixes).

### Default Input Distributions

Instances consist of both synthetic and real strings.

The large instances are:

- A trigram string of length 100 Million, generated with:  
`trigramString 100000000 <filename>`

- `chr22.dna` is a DNA sequence.  It consists only of the
characters C,G,C,A,N and has about 34 million characters.

- `etext99` is text from the project Guttenberg.  It has about 105 Million characters.

- `wikisamp.xml` is a sample from wikipedia's xml source files.  It has
exactly 100 million characters.

The small instances are:

- A trigram string of length 10 Million, generated with:  
`trigramString 10000000 <filename>`

- `chr22.dna` as for the large instances.

### Input and Output File Formats

The input needs to be a file of characters (no null characters).
The output needs to be in the [sequence file
format](../fileFormats/sequence.html) with integer type.

