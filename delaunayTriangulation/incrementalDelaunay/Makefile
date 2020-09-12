include parallelDefs

# required files
NEAREST_NEIGHBORS = geometry.h nearestNeighbors.h blockRadixSort.h quickSort.h transpose.h
DELAUNAY = topology.h delaunayQs.h delaunayDefs.h delaunay.C 
GLOBAL_REQUIRE = parallel.h sequence.h $(NEAREST_NEIGHBORS) $(DELAUNAY)
LOCAL_REQUIRE = octTree.h
OBJS = delaunay.o 

include ../common/timeRequiredFiles
include MakeBench
