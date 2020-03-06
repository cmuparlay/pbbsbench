# pbbsbench
New version of pbbs benchmarks (work in progress)

This repository uses a submodule (pbbslib).    To pull the submodule:

\> git submodule init

<prompt> git submodule update

To run a test, for example, try:

<prompt> cd comparisonSort/sampleSort

<prompt> make

<prompt> ./testInputs -x -r 3

The -x means don't check correctness, the -r 3 means do three trials.
