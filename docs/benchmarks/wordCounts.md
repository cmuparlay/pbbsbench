---
title: Word Counts
---

# Word Counts (WC)

Given an input string counts the number of occurences of each word in
the string.  All non characters (not in [a-z] or [A-Z]) are replaced
with a blank, and all upper-case characters are converted to lower
case.  A word is a maximal contiguous sequence of characters in the
resulting string.

The output is a sequence of pairs, each consisting of a string (word)
and a count of how many times it appears. Ordering does not matter.

### Default Input Distributions

Instances consist of both synthetic and real strings.

The large instances are:

- A trigram string of length 250 Million, generated with:  
`trigramString 250000000 <filename>`

- `etext99` is text from the project Guttenberg.  It has about 105
Million characters.

- `wikipedia250M.txt` is a sample from wikipedia's xml source files.  It has 
exactly 250 million characters.

The small instances are:

- A trigram string of length 25 Million, generated with:  
`trigramString 25000000 <filename>`

- `wikisamp.xml` is a sample from wikipedia's xml source files.  It has 
exactly 100 million characters.

### Input and Output File Formats

The input is a text file and output need to be in the [sequence file format](../fileFormats/sequence.html),
with type `StringIntPair`.

