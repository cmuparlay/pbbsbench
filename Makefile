
DEFAULT_BENCHMARKS = integerSort/parallelRadixSort comparisonSort/sampleSort comparisonSort/serialSort removeDuplicates/serial_hash removeDuplicates/parlayhash histogram/parallel histogram/sequential wordCounts/histogram wordCounts/serial invertedIndex/sequential invertedIndex/parallel suffixArray/parallelRange suffixArray/serialDivsufsort longestRepeatedSubstring/doubling classify/decisionTree minSpanningForest/parallelFilterKruskal minSpanningForest/serialMST spanningForest/ndST spanningForest/serialST breadthFirstSearch/backForwardBFS breadthFirstSearch/serialBFS maximalMatching/serialMatching maximalMatching/incrementalMatching maximalIndependentSet/ndMIS maximalIndependentSet/serialMIS nearestNeighbors/octTree rayCast/kdTree convexHull/quickHull convexHull/serialHull delaunayTriangulation/incrementalDelaunay delaunayRefine/incrementalRefine rangeQuery2d/parallelPlaneSweep rangeQuery2d/serial nBody/parallelCK

EXT_BENCHMARKS = comparisonSort/quickSort comparisonSort/mergeSort comparisonSort/stableSampleSort comparisonSort/ips4o removeDuplicates/serial_sort suffixArray/parallelKS spanningForest/incrementalST breadthFirstSearch/simpleBFS breadthFirstSearch/deterministicBFS maximalIndependentSet/incrementalMIS 

DATA_GENERATORS = sequenceData graphData geometryData

all : $(DEFAULT_BENCHMARKS) $(DATA_GENERATORS)

ext : all $(EXT_BENCHMARKS)

$(DEFAULT_BENCHMARKS) : FORCE
	cd $@; make -s

$(EXT_BENCHMARKS) : FORCE
	cd $@; make -s

$(DATA_GENERATORS) : FORCE
	cd testData/$@; make -j -s

FORCE :
