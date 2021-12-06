---
title: Sequence File Formats
---

# Sequence File Format

The sequence file format is used for representing a sequence of values.
The files are stored in ascii and need to be in the format:

```
sequence<type>
<v0>
<v1>
...
<v_(n-1)>
```
Each item after the `sequence<type>`
declaration represents one element of the sequence and
`<type>` specifies the type of each elements.
The items are delimited by any consecutive sequence of
delimiter characters: **tab**, **space**,
**line feed** (ascii 0x0A), and **carriage
return** (ascii 0x0D).  The `<type>` can
be one of the following:

- `Int`: integers in decimal notation.  A minus sign (-) is used to indicate a negative integer.   
- `Double`: double precision floating point numbers in decimal or exponential notation.
- `String`:  character strings (can contain any sequence
of non-delimiter character).   There is no limit on the length of the strings.
- `<type1><type2>Pair`:
where `<type1>` and `<type2>` are any
of the types above.  This represents a pair of values of the given
types separated by a delimiter character.  For example:
    > `sequenceStringIntPair`  
    > `this 22`  
    > `is    444`  
    > `a 12345`  
    > `sequence 11`  
    > `of -22 nine`  
    > `14`  
    > `string 1`  
    >
    > `integer 16`  
    > `pairs 0`  

The missing and extra carriage returns are meant to illustrate that
there is no distinction between the delimiting characters.

Files can start and end with delimiters, which are ignored.
