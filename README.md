# pbbsbench
New version of pbbs benchmarks (work in progress)

This repository uses a submodule (pbbslib).    To pull the submodule:

\> git submodule init

\> git submodule update

To run a test, for example, try:

\> cd comparisonSort/sampleSort

\> make

\> ./testInputs -x -r 3

The -x means don't check correctness, the -r 3 means do three trials.
