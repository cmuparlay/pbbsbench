include common/parallelDefs

COMMON = common/graph.h common/graphIO.h common/graphUtils.h
GENERATORS = rMatGraph gridGraph randLocalGraph nBy2Comps lineGraph addWeights adjToEdgeArray edgeArrayToAdj 

NOTUPDATED_GENERATORS = powerGraph addWeights randDoubleVector fromAdjIdx adjElimSelfEdges starGraph combGraph adjGraphAddWeights binTree randGraph reorderGraph randomizeGraphOrder adjGraphAddSourceSink dimacsToFlowGraph adjToBinary adjWghToBinary

.PHONY: all clean
all: $(GENERATORS)

%.o : %.C $(COMMON)
	$(CC) $(CFLAGS) -c $< -o $@

rMatGraph : rMatGraph.o 
	$(CC) $(LFLAGS) -o $@ rMatGraph.o 

gridGraph : gridGraph.o 
	$(CC) $(LFLAGS) -o $@ gridGraph.o 

nBy2Comps : nBy2Comps.o 
	$(CC) $(LFLAGS) -o $@ nBy2Comps.o 

lineGraph : lineGraph.o 
	$(CC) $(LFLAGS) -o $@ lineGraph.o 

powerGraph : powerGraph.o 
	$(CC) $(LFLAGS) -o $@ powerGraph.o 

randLocalGraph : randLocalGraph.o 
	$(CC) $(LFLAGS) -o $@ randLocalGraph.o 

randGraph : randGraph.o 
	$(CC) $(LFLAGS) -o $@ randGraph.o 

addWeights : addWeights.o 
	$(CC) $(LFLAGS) -o $@ addWeights.o

fromAdjIdx : fromAdjIdx.o 
	$(CC) $(LFLAGS) -o $@ fromAdjIdx.o

randDoubleVector : randDoubleVector.o
	$(CC) $(LFLAGS) -o $@ randDoubleVector.o

adjToEdgeArray : adjToEdgeArray.C $(COMMON)
	$(CC) $(CFLAGS) $(LFLAGS) -o $@ adjToEdgeArray.C

adjElimSelfEdges : adjElimSelfEdges.C $(COMMON)
	$(CC) $(CFLAGS) $(LFLAGS) -o $@ adjElimSelfEdges.C

edgeArrayToAdj : edgeArrayToAdj.C $(COMMON)
	$(CC) $(CFLAGS) $(LFLAGS) -o $@ edgeArrayToAdj.C

starGraph : starGraph.o 
	$(CC) $(LFLAGS) -o $@ starGraph.o 

combGraph : combGraph.o 
	$(CC) $(LFLAGS) -o $@ combGraph.o 

expGraph : expGraph.o 
	$(CC) $(LFLAGS) -o $@ expGraph.o 

binTree : binTree.o 
	$(CC) $(LFLAGS) -o $@ binTree.o 

reorderGraph : reorderGraph.o 
	$(CC) $(LFLAGS) -o $@ reorderGraph.o 

randomizeGraphOrder : randomizeGraphOrder.o 
	$(CC) $(LFLAGS) -o $@ randomizeGraphOrder.o 

adjGraphAddWeights : adjGraphAddWeights.o 
	$(CC) $(LFLAGS) -o $@ adjGraphAddWeights.o

adjGraphAddSourceSink : adjGraphAddSourceSink.o
	$(CC) $(LFLAGS) -o $@ adjGraphAddSourceSink.o

dimacsToFlowGraph : dimacsToFlowGraph.o
	$(CC) $(LFLAGS) -o $@ dimacsToFlowGraph.o

flowGraphToDimacs : flowGraphToDimacs.o
	$(CC) $(LFLAGS) -o $@ flowGraphToDimacs.o

adjToBinary : adjToBinary.o 
	$(CC) $(LFLAGS) -o $@ adjToBinary.o

adjWghToBinary : adjWghToBinary.o 
	$(CC) $(LFLAGS) -o $@ adjWghToBinary.o

clean :
	rm -f *.o $(GENERATORS)
#causes error
#   	cd maxFlowGens; make clean
#	make clean -s -C data

cleansrc : 
	make -s clean
	rm -f $(COMMON) 
