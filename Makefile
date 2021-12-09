
DEFAULT_BENCHMARKS = integerSort/parallelRadixSort comparisonSort/sampleSort comparisonSort/serialSort removeDuplicates/serial_hash removeDuplicates/parlayhash histogram/parallel histogram/sequential wordCounts/histogram wordCounts/serial invertedIndex/sequential invertedIndex/parallel suffixArray/parallelRange suffixArray/serialDivsufsort longestRepeatedSubstring/doubling classify/decisionTree minSpanningForest/parallelFilterKruskal minSpanningForest/serialMST spanningForest/ndST spanningForest/serialST breadthFirstSearch/backForwardBFS breadthFirstSearch/serialBFS maximalMatching/serialMatching maximalMatching/incrementalMatching maximalIndependentSet/ndMIS maximalIndependentSet/serialMIS nearestNeighbors/octTree rayCast/kdTree convexHull/quickHull convexHull/serialHull delaunayTriangulation/incrementalDelaunay delaunayRefine/incrementalRefine rangeQuery2d/parallelPlaneSweep rangeQuery2d/serial nBody/parallelCK

DATA_GENERATORS = sequenceData graphData geometryData

all : $(DEFAULT_BENCHMARKS) $(DATA_GENERATORS)

$(DEFAULT_BENCHMARKS) : FORCE
	cd $@; make -s

$(DATA_GENERATORS) : FORCE
	cd testData/$@; make -j -s

FORCE :
