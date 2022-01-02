---
title: The PBBS Benchmark Suite (V2)
---

#  The PBBS Benchmark Suite (V2)

This documents the Problem Based Benchmark Suite (PBBS), a collection
of over 20 benchmarks defined in terms of their IO characteristics.
They are designed to make it possible to compare different algorithms,
or implementations in different programming languages.

The current benchmarks are [**listed here**](benchmarks/index.html).

For each benchmark the suite provides:

- the specification of the input and expected output for the problem,
- the specification of a default set of input instances,
- code for generating inputs (written to a file),
- code for checking correctness of output (read from a file),
- code for timing the benchmark across the instances,
- a default parallel implementation,
- a default sequential implementation (for most benchmarks),
- a variety of other implementations (for some benchmarks).

The benchmarks are designed to be agnostic to the programming
language.  However, the framework is mostly in C++ and some of the
tools are easier to use with C++.  For example there is a timing
driver for C++ that can be linked with any implementation of a
benchmark.

### Table of Contents

1. [Getting Started](#getting-started)
   - [Requirements](#requirements)
   - [Running the Benchmarks](#running-the-benchmarks)
2. [More Details](#more-details)
   - [Options for runall](#options-for-runall)
   - [The Benchmark Directories](#the-benchmark-directories)
   - [Input Instances and Data Generators](#input-instances-and-data-generators)
   - [Timing the Benchmarks](#timing-the-benchmarks)
   - [The Driver](#the-driver)
   - [Checking Correctness](#checking-correctness)
   
## Getting Started

The benchmark suite can be downloaded from github using:

```
> git clone https://github.com/cmuparlay/pbbsbench.git
```

It uses two submodules, which can be initialized with:

```
> git submodule update --init
```

### Requirements

The benchmarks have been tested on Ubuntu and MacOS.

The software requirements are:

- C++-17 compiler (tested with gcc and clang)

The system requirements are

- for small data (12GB of RAM and 2GB of disk)
- for large data (64GB of RAM and 10GB of disk)

The following are not required, but will give better performance

- `jemalloc`  (only gives slight performance improvement)
- 20+ cores (the more cores the faster)
- numactl installed (if this is not installed you need to run `./runall -nonuma`)

### Running the benchmarks

The command `./runall` will compile and run all the benchmarks
reported but will take a couple hours.  For a faster run, try:

```
  ./runall -par -small
```
  
This only compiles and runs the parallel benchmarks, which run much faster, and on
significantly smaller input data (an order of magnitude smaller for some benchmarks).
This runs in 12 minutes on a 20 core (40 hyperthread) machine.

You can also test individual benchmarks.   For example, you can test the
parallel comparison sort using:

```
  ./runall -only comparisonSort/sampleSort
 ```
  
This will run the parallel sampleSort on the default (full sized) inputs.
The call

```
  ./runall -small -only comparisonSort/sampleSort
```
  
will run it on the smaller inputs.  More details on arguments are
given below.

## More Details

### Options for `runall`

The ./runall has the following options which can be extracted by using
`./runall -h`.

```
  -scale    : this runs it on a range of different thread counts up the the number of threads on the machine
  -small    : runs tests on smaller inputs (calls ./testInput_small instead of ./testInput).
  -par      : only run benchmarks that are parallel (saves time)
  -only <name>   : only run a particular benchmark
  -notime   : only compile the benchmarks
  -nonuma   : don't use numactl
  -nocheck  : don't check correctness of results (saves time)
```
  
For the `-only` option use the path to the implementation, e.g.

```
  ./runall -only comparisonSort/sampleSort
```

### The Benchmarks Directories

Within the `benchmarks` directory at toplevel is a subdirectory
for each of the benchmarks.
The [benchmarks page](benchmarks/index.html) gives a listing of these
the benchmarks along with their specification.

Within each benchmarks is a subdirectory for each of the
implementations.  Each benchmark also has some directories shared
across implementations.  In particular each has a directory called
`bench` which contains the driver code, testing code, and
specification of default inputs.  Each benchmark also has a xxxData
page containing data generators for the benchmark (here xxx varies by
benchmark).

Within each implementation directory, you can run `make` to make the
executable, and then run `./testInputs` to run the benchmarks.  These
are run automatically by the ./runall script.  On a machine with
multiple chips, using `numactl -i all ./testInputs` will give better
results.  `./testInputs_small` will use the smaller inputs.

The `testInputs` script has several options including:

```
  -x : do not check the output
  -r <count>  : number of rounds to use
  -p <count>  : number of threads to use
  ```
  
The actual inputs are specified in the script and can be changed if desired.

### Input Instances and Data Generators

Each benchmark has suggested input instances.   There are two sets of
such instances, one small, and the other large.   Many instances are
synthetic and generated by one of our data generators.  All the data
generators are in the toplevel directory `testData`.   There are three
subdirectories `sequenceData`, `graphData` and `geometryData`, each
with its own generators.   The various testing scripts (e.g., `runAll`
and `testInputs`) compile and run the generators so there should be no
need to generate the data by hand.    By default 
generated test inputs are deleted after reuse.

### Timing the Benchmarks

Users can use the benchmark suite as they please, but here are the
expectations for timing an implementation so that it can be properly
compared with other implementations, including the default ones.  We
do not expect a full end-to-end timing of the executable because often
the time is dominated by reading the input file and possibly writing
an output file, which have little to do with the benchmark.  Instead
we expect a benchmark to wrap a timer (real time) around the benchmark
itself.

The reading of the data should be before the timing starts, and at the
start the input data should be of the type suggested by the benchmark
--- e.g. a sequence, or graph.  In general it does not matter how this
input is formatted.  For graphs, for example, either adjacency-arrays
or edge arrays are fine.  However, the data cannot be pre-processed in
any significant way.  For example, the edges of a weighted graph
cannot be sorted by weight (this would make MST much faster), and the
vertices should not be reordered for locality.  If someone is
interested in studying the benefit of reordering, for example, we
suggest they create new input instances that are reordered how they
want.  In this way they could be compared with other code in a fair
way.

It is considered fine to run the benchmark a few times before starting
the timing.  This can ensure, for example, that virtual memory is
paged in.  However, it is of course important the runs do no leave
problem specific information (e.g. cached partial results) that would
benefit future runs.

The result type of the benchmark should be as described in the
benchmark specification--e.g., a sequence, or a string.  As with the
input there is flexibility on representation.  However, the result
cannot be "lazy" in that it is not fully computed, and to be completed
when output to a file.  Outputting the result to a file for the purpose
of checking it, or other checking of correctness within the benchmark,
should be outside of the timing.

### The Driver

A benchmark implementation should generate a driver executable.
The driver executable should allow the user to specify an input file,
an optional output file, and how many times they want to run it.
To use it with our `testInputs` script the command line should be of
the form:

```
  <bench> [-o <outfile>] [-r <numrounds>] <infile>
```

Other arguments are welcome.
Each benchmark expects the driver executable `<bench>` to 
have a specific name, e.g. `sort` for comparison sort.

To use with the testInput script the driver should report times to
stdout, one per run on its own line in the format:

`PBBS Time: <seconds>`

Any other output is ignored. We supply C++ driver code so if written
in C++ the user can just define a function with the specified input
format and output format.  If they want to use a different input
format, they will have to modify the driver.

### Checking Correctness

Most benchmarks come with programs that test for correctness.  Some
where there is no absolute correctness criteria, report some measure
of acccuracy instead.  For example the `classify` checker reports the
percentage of the test vectors that were correctly classified.  The
checking is file based (to allow cross language communication).  Most
checkers take the original input file, as well as the output file
generated by the benchmark, and check the two against each other,
reporting any errors it finds. If nothing is reported, the test passed.
We cannot guarantee that all checkers completely checks for
correctness.
