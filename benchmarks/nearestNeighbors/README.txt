INTRO
This folder contains a number of parallel nearest neighbor implementations. Each is contained in its own folder. To run it, navigate to the folder and run 
the "make" command. Then, calling "./testInputs" runs some tests on various data files. To run a nearest neighbor algorithm on an arbitrary data file in D
dimensions and searching for K nearest neighbors, execute "./neighbors -d D -k K /path/to/data." Most nearest neighbor algorithms have a "report_stats" 
parameter in the "neighbors.h" file that can be turned on or off. 

Currently only 2 and 3 dimensions are supported. Calculating up to 100 nearest neighbors is supported.

OCTTREE
This is our homemade implementation. It uses a kd-tree where the splitting rule is based on the Morton ordering of the point set. It has three options for 
traversing the kd-tree (this is the parameter "algorithm_version"). The "root-based" version starts nearest neighbor searches from the root of the kd-tree, 
and is applicable when searching points already in the kd-tree, or points not in the kd-tree, which we refer to as dynamic queries. The "bit-based" version 
navigates to the leaf of the kd-tree and then begins the nearest neighbor search starting from the leaf; it works for dynamic or non-dynamic queries. The 
"map-based" version only works for non-dynamic queries; it stores pointers from each point to the leaf it is contained in, and begins the nearest neighbor
search from the leaf.  

The kd-tree also supports batch-dynamic insertions and deletions. 

All these features are demonstrated in the "neighbors.h" file in the OCTREE folder.

-------BENCHMARKS---------

NAIVE
This is a very basic O(n^2) solution. Only works for k=1.

STANN
This implementation is based on the paper "Fast construction of k-nearest neighbor graphs for point clouds" by Michael Connor and Piyush Kumar. Most of the
code is from those authors. There are four different benchmarks within STANN. Two are specialized for computing the k-nearest neighbor graph of a point
set. Of those, "KNNG" is Connor and Kumar's implementation, which uses parallel primitives from OpenMP. "ParlayKNNG" is a version we modified to use
ParlayLib's parallel primitives, as well as our own sorting algorithms, and some changes to the base case. The other two are specialized for dynamic 
queries. "KNN" and "ParlayKNN" are Connor and Kumar's implementation, and our own, respectively.

CHAN05
This implementation is a parallel version of the code provided in Timothy Chan's paper "A minimalist's implementation of an approximate nearest neighbor 
algorithm in fixed dimensions." It only works for k=1.

CKNN
This is a parallel version of the algorithm found in Paul B. Callahan and S. Rao Koseraju's paper "A Decomposition of Multidimensional Point Sets with
Applications to k-Nearest-Neighbors and n-Body Potential Fields." It works for arbitrary k, but is not very practical for k>30.

CGAL
This is the kd-tree based nearest neighbor algorithm from the Computational Geometry Algorithms Library. It is under construction with numerous issues. 
The code itself has some contention, and due to this we recommend running only on one chip. To run it, the user needs to download CGAL and change 
the compiler flags to point to it correctly. To compile, it also currently requires the "#ifndef" in line 236 of "alloc.h" to be changed to 
"#ifdef".

