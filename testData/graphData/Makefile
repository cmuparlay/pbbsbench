include common/parallelDefs

COMMON = common/graph.h common/graphIO.h common/graphUtils.h
GENERATORS = rMatGraph gridGraph randLocalGraph nBy2Comps lineGraph addWeights adjToEdgeArray edgeArrayToAdj

NOTUPDATED_GENERATORS = powerGraph addWeights randDoubleVector fromAdjIdx adjElimSelfEdges starGraph combGraph adjGraphAddWeights binTree randGraph reorderGraph randomizeGraphOrder adjGraphAddSourceSink dimacsToFlowGraph adjToBinary adjWghToBinary

.PHONY: all clean
all: $(GENERATORS)

%.o : %.C $(COMMON)
	$(CC) $(CFLAGS) -c $< -o $@

rMatGraph : rMatGraph.o
	$(CC) -o $@ rMatGraph.o  $(LFLAGS)

gridGraph : gridGraph.o
	$(CC) -o $@ gridGraph.o  $(LFLAGS)

nBy2Comps : nBy2Comps.o
	$(CC) -o $@ nBy2Comps.o  $(LFLAGS)

lineGraph : lineGraph.o
	$(CC) -o $@ lineGraph.o  $(LFLAGS)

powerGraph : powerGraph.o
	$(CC) -o $@ powerGraph.o  $(LFLAGS)

randLocalGraph : randLocalGraph.o
	$(CC) -o $@ randLocalGraph.o  $(LFLAGS)

randGraph : randGraph.o
	$(CC) -o $@ randGraph.o  $(LFLAGS)

addWeights : addWeights.o
	$(CC) -o $@ addWeights.o $(LFLAGS)

fromAdjIdx : fromAdjIdx.o
	$(CC) -o $@ fromAdjIdx.o $(LFLAGS)

randDoubleVector : randDoubleVector.o
	$(CC) -o $@ randDoubleVector.o $(LFLAGS)

adjToEdgeArray : adjToEdgeArray.C $(COMMON)
	$(CC) $(CFLAGS) -o $@ adjToEdgeArray.C $(LFLAGS)

adjElimSelfEdges : adjElimSelfEdges.C $(COMMON)
	$(CC) $(CFLAGS) -o $@ adjElimSelfEdges.C $(LFLAGS)

edgeArrayToAdj : edgeArrayToAdj.C $(COMMON)
	$(CC) $(CFLAGS) -o $@ edgeArrayToAdj.C $(LFLAGS)

starGraph : starGraph.o
	$(CC) -o $@ starGraph.o  $(LFLAGS)

combGraph : combGraph.o
	$(CC) -o $@ combGraph.o  $(LFLAGS)

expGraph : expGraph.o
	$(CC) -o $@ expGraph.o  $(LFLAGS)

binTree : binTree.o
	$(CC) -o $@ binTree.o  $(LFLAGS)

reorderGraph : reorderGraph.o
	$(CC) -o $@ reorderGraph.o  $(LFLAGS)

randomizeGraphOrder : randomizeGraphOrder.o
	$(CC) -o $@ randomizeGraphOrder.o  $(LFLAGS)

adjGraphAddWeights : adjGraphAddWeights.o
	$(CC) -o $@ adjGraphAddWeights.o $(LFLAGS)

adjGraphAddSourceSink : adjGraphAddSourceSink.o
	$(CC) -o $@ adjGraphAddSourceSink.o $(LFLAGS)

dimacsToFlowGraph : dimacsToFlowGraph.o
	$(CC) -o $@ dimacsToFlowGraph.o $(LFLAGS)

flowGraphToDimacs : flowGraphToDimacs.o
	$(CC) -o $@ flowGraphToDimacs.o $(LFLAGS)

adjToBinary : adjToBinary.o
	$(CC) -o $@ adjToBinary.o $(LFLAGS)

adjWghToBinary : adjWghToBinary.o
	$(CC) -o $@ adjWghToBinary.o $(LFLAGS)

clean :
	rm -f *.o $(GENERATORS)
#causes error
#   	cd maxFlowGens; make clean
#	make clean -s -C data

cleansrc :
	make -s clean
	rm -f $(COMMON)
