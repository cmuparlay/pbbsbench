include parallelDefs

COMMON = dataGen.h utils.h IO.h parseCommandLine.h graph.h graphIO.h graphUtils.h parallel.h sequence.h blockRadixSort.h deterministicHash.h transpose.h quickSort.h
GENERATORS = rMatGraph gridGraph nBy2Comps lineGraph powerGraph randLocalGraph addWeights randDoubleVector fromAdjIdx adjToEdgeArray adjElimSelfEdges edgeArrayToAdj starGraph combGraph adjGraphAddWeights binTree randGraph reorderGraph randomizeGraphOrder adjGraphAddSourceSink dimacsToFlowGraph adjToBinary adjWghToBinary

.PHONY: all clean
all: $(GENERATORS)

$(COMMON) :
	ln -s ../../common/$@ .

%.o : %.C $(COMMON)
	$(PCC) $(PCFLAGS) -c $< -o $@

rMatGraph : rMatGraph.o 
	$(PCC) $(PLFLAGS) -o $@ rMatGraph.o 

gridGraph : gridGraph.o 
	$(PCC) $(PLFLAGS) -o $@ gridGraph.o 

nBy2Comps : nBy2Comps.o 
	$(PCC) $(PLFLAGS) -o $@ nBy2Comps.o 

lineGraph : lineGraph.o 
	$(PCC) $(PLFLAGS) -o $@ lineGraph.o 

powerGraph : powerGraph.o 
	$(PCC) $(PLFLAGS) -o $@ powerGraph.o 

randLocalGraph : randLocalGraph.o 
	$(PCC) $(PLFLAGS) -o $@ randLocalGraph.o 

randGraph : randGraph.o 
	$(PCC) $(PLFLAGS) -o $@ randGraph.o 

addWeights : addWeights.o 
	$(PCC) $(PLFLAGS) -o $@ addWeights.o

fromAdjIdx : fromAdjIdx.o 
	$(PCC) $(PLFLAGS) -o $@ fromAdjIdx.o

randDoubleVector : randDoubleVector.o
	$(PCC) $(PLFLAGS) -o $@ randDoubleVector.o

adjToEdgeArray : adjToEdgeArray.C $(COMMON)
	$(PCC) $(PCFLAGS) $(PLFLAGS) -o $@ adjToEdgeArray.C

adjElimSelfEdges : adjElimSelfEdges.C $(COMMON)
	$(PCC) $(PCFLAGS) $(PLFLAGS) -o $@ adjElimSelfEdges.C

edgeArrayToAdj : edgeArrayToAdj.C $(COMMON)
	$(PCC) $(PCFLAGS) $(PLFLAGS) -o $@ edgeArrayToAdj.C

starGraph : starGraph.o 
	$(PCC) $(PLFLAGS) -o $@ starGraph.o 

combGraph : combGraph.o 
	$(PCC) $(PLFLAGS) -o $@ combGraph.o 

expGraph : expGraph.o 
	$(PCC) $(PLFLAGS) -o $@ expGraph.o 

binTree : binTree.o 
	$(PCC) $(PLFLAGS) -o $@ binTree.o 

reorderGraph : reorderGraph.o 
	$(PCC) $(PLFLAGS) -o $@ reorderGraph.o 

randomizeGraphOrder : randomizeGraphOrder.o 
	$(PCC) $(PLFLAGS) -o $@ randomizeGraphOrder.o 

adjGraphAddWeights : adjGraphAddWeights.o 
	$(PCC) $(PLFLAGS) -o $@ adjGraphAddWeights.o

adjGraphAddSourceSink : adjGraphAddSourceSink.o
	$(PCC) $(PLFLAGS) -o $@ adjGraphAddSourceSink.o

dimacsToFlowGraph : dimacsToFlowGraph.o
	$(PCC) $(PLFLAGS) -o $@ dimacsToFlowGraph.o

flowGraphToDimacs : flowGraphToDimacs.o
	$(PCC) $(PLFLAGS) -o $@ flowGraphToDimacs.o

adjToBinary : adjToBinary.o 
	$(PCC) $(PLFLAGS) -o $@ adjToBinary.o

adjWghToBinary : adjWghToBinary.o 
	$(PCC) $(PLFLAGS) -o $@ adjWghToBinary.o

clean :
	rm -f *.o $(GENERATORS)
#causes error
#   	cd maxFlowGens; make clean
#	make clean -s -C data

cleansrc : 
	make -s clean
	rm -f $(COMMON) 
