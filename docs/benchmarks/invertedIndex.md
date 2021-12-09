---
title: Inverted Index
---

# Inverted Index (IIDX)

Generates an inverted index from a text file.  It assumes each
document starts with a given "start string".  In this benchmark the
start string is `<doc` (as used in wikipedia entries).  The benchmark
breaks the file into "documents" (strings) based on this document
start string.  The start strings are dropped, and documents are given
integer identifiers starting from 0 and in the order they appear in
the file.

For each document all upper case letters (A-Z) are converted to
lowercase all lowercase letters are kept and all others are blacked
out.  What remains is broken into words (i.e. contiguous sequences of
letters).  We now build an "inverted index", where each entry consists
of a word and all the documents it belongs in.  The entries are are
sorted lexicographically by word, and the list of entries are sorted
by document identifier.  The output is then a string containing one
entry per line, starting with the word, and then the document
identifiers (in ascii).

For example, the input string

    <doc this is a+String <doc and_a2 sTrIng

 would generate the inverted index string:

```
    a 0 1
    and 1
    is 0
    this 0
    string 0 1
  ```
  
The timing must include the full cost of from parsing the input string to generating the final output string.


### Default Input Distributions

The test distributions are two files taken from wikipedia.

- `wikisamp.xml` : these are just summaries of each entry and short documents

- `wikipedia250M.txt` : these are full documents

### Input and Output File Formats 

The input is an ascii string containing the documents.   The output is
as ascii string as described above. 
