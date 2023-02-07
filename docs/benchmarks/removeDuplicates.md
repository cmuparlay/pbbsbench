---
title: Remove Duplicates
---

# Remove Duplicates (DDUP)

Given a sequence of elements of any (uniform) type, remove any
duplicates returning only one element for each key value. Each element
consists of a key and possibly auxiliary data and the implementation
must be based on the following three user supplied functions:

* A hash function that maps each key to an unsigned integer,
* a comparison function on the keys defining a total order, and
* a comparison function on the auxiliary data defining a partial order.

The comparison function on the auxiliary data is used to decide which
element to keep when multiple elements have the same key--a maximal
value with respect to the partial order must be kept. When there is no
auxiliary data this function should always returns false. The code
must not take advantage of the specific key and auxiliary types beyond
the hash and comparison function.

### Default Input Distributions

The distributions are as follows:

  * Unsigned integers generated uniformly at random in the range [0:n).

        randomSeq -t int <n> <filename>

  * Unsigned integers generated at random with repeats appearing in an
    exponential distribution.

        exptSeq -t int <n> <filename>

  * Strings from a tri-gram distribution.

        trigramSeq <n> <filename>

### Input and Output File Formats

The input and output should be in the [sequence file
format](../fileFormats/sequence.html) both with the same element
types. It output must contain one element for each key in the input,
and if there is auxiliary data, a maximal auxiliary value for that
key. The output can be ordered in any way.
