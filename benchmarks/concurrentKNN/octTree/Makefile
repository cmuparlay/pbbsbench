include common/parallelDefs

BENCH = neighbors
REQUIRE = oct_tree.h k_nearest_neighbors.h
CFLAGS += -DNoHelp #-DVersioned -DLazyStamp #-DHandOverHand

include common/MakeBench

lucy_utils: lucy_utils.cpp
	$(CC) $(CFLAGS) -o lucy_utils lucy_utils.cpp $(LFLAGS) 