---
title: N-body Force Calculation
---

# N-body Force Calculation (NBODY)

Given n points in 3 dimensions calculate the gravitational force
vector on each point due to all other points.  This force vector can
be an approximation.  For this benchmark we use a gravitational
constant of 1 giving the force between two particles located
at **r**<sub>1</sub> and **r**<sub>2</sub> as:

> <strong>F</strong><sub>12</sub> =
m<sub>1</sub> m<sub>2</sub> (<strong>r</strong><sub>2</sub> -
<strong>r</strong><sub>1</sub>)/ \|\|(<strong>r</strong><sub>2</sub> -
<strong>r</strong><sub>1</sub>)\|\|<sup>3</sup>.

The input is a sequence of 3-d points, each coordinate which is a
double-precision floating-point number.    The output is a sequence of
equal length containing a 3-d vector for each of the points, again in
double-precision floating-point precision for each coordinate.
We say a reported force F' has error bound epsilon if

> \|\| F - F' \|\| / \|\| F \|\| < epsilon,

where F is the actual force.  By default all reported forces must have
an error bound of epsilon = 10<sup>-6</sup>.  Ideally there should be
a input parameter that gives an accuracy, time tradeoff.

### Default Input Distributions

The default distributions are the following:

- Points chosen uniformly at random on the surface of a unit sphere.   Should be 
generated with:  
`randPoints -S -d 3 <n> <filename>`

- Points chosen uniformly at random within a unit cube.   Should be 
generated with:  
`randPoints  -d 3 <n> <filename>`

- Points chosen at random from the Plummer distribution.   Should be 
generated with:  
`randPoints -k -d 3 <n> <filename>`

The large size is n = 1 million, and the small size is n = .1 million.

### Input and Output File Formats

The input and output need to be in the [3dpoints file format](../fileFormats/geometry.html#points).
