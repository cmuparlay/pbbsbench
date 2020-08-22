// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#pragma once

#include "parlay/parallel.h"
#include "common/geometry.h"
#include "common/dataGen.h"

template <class coord>
point2d<coord> rand2d(size_t i) {
  size_t s1 = i;
  size_t s2 = i + dataGen::hash<size_t>(s1);
  return point2d<coord>(2*dataGen::hash<double>(s1)-1,
			2*dataGen::hash<double>(s2)-1);
}

template <class coord>
point3d<coord> rand3d(size_t i) {
  size_t s1 = i;
  size_t s2 = i + dataGen::hash<size_t>(s1);
  size_t s3 = 2*i + dataGen::hash<size_t>(s2);
  return point3d<coord>(2*dataGen::hash<double>(s1)-1,
			2*dataGen::hash<double>(s2)-1,
			2*dataGen::hash<double>(s3)-1);
}

template <class coord>
point2d<coord> randInUnitSphere2d(size_t i) {
  size_t j = 0;
  vector2d<coord> v;
  do {
    size_t o = dataGen::hash<size_t>(j++);
    v = vector2d<coord>(rand2d<coord>(o+i));
  } while (v.Length() > 1.0);
  return point2d<coord>(v);
}

template <class coord>
point3d<coord> randInUnitSphere3d(size_t i) {
  size_t j = 0;
  vector3d<coord> v;
  do {
    size_t o = dataGen::hash<size_t>(j++);
    v = vector3d<coord>(rand3d<coord>(o+i));
  } while (v.Length() > 1.0);
  return point3d<coord>(v);
}

template <class coord>
point2d<coord> randOnUnitSphere2d(size_t i) {
  vector2d<coord> v = vector2d<coord>(randInUnitSphere2d<coord>(i));
  return point2d<coord>(v/v.Length());
}

template <class coord>
point3d<coord> randOnUnitSphere3d(size_t i) {
  vector3d<coord> v = vector3d<coord>(randInUnitSphere3d<coord>(i));
  return point3d<coord>(v/v.Length());
}
