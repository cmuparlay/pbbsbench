---
title: Longest Repeated Substring
---

# Longest Repeated Substring (LRS)

Given a string identifies the longest repeated substrings.   The two
repeated substrings may overlap.   The input is a string with no null
characters (but any other value is OK). The output is just a tripple
consisting of the length of the string, the start position of one
occurence, and the start position of another occurence.   Start
positions are zero based.

### Default Input Distributions

Instances consist of both synthetic and real strings.
The large instances are:

- `chr22.dna` is a DNA sequence.  It consists only of the
characters C,G,C,A,N and has about 34 million characters.

- `etext99` is text from the project Guttenberg.  It has about 105
Million characters.

- `wikisamp.xml` is a sample from wikipedia's xml source files.  It has
exactly 100 million characters.

The small instances only includes `chr22.dna`

### Input and Output File Formats

The input needs to be a file of characters (no null characters).
The output is simply an ascii file with three numbers in it: the
length, and the two positions.

