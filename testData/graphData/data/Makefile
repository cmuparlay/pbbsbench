GENERATORS = ../gridGraph ../rMatGraph ../randLocalGraph ../addWeights ../randDoubleVector ../starGraph ../combGraph ../expGraph ../adjGraphAddWeights ../adjGraphAddSourceSink ../dimacsToFlowGraph ../flowGraphToDimacs

$(GENERATORS) :
	cd ..; make -s $(subst ../,,$@)

MAXFLOWGENS = ../maxFlowGens/bin/ac ../maxFlowGens/bin/ak ../maxFlowGens/bin/genrmf ../maxFlowGens/bin/washington ../maxFlowGens/bin/genrmf-wrapper ../maxFlowGens/bin/washington-wrapper
$(MAXFLOWGENS):
	cd ../maxFlowGens; make

clean :
	rm -f *0*

2Dgrid_E_% : ../gridGraph
	../gridGraph -d 2 $(subst 2Dgrid_E_,,$@) $@

3Dgrid_E_% : ../gridGraph
	../gridGraph -d 3 $(subst 3Dgrid_E_,,$@) $@

rMatGraph_E_5_% : ../rMatGraph
	../rMatGraph  -a .5 -b .1 -m `expr 5 \* $(subst rMatGraph_E_5_,,$@)` $(subst rMatGraph_E_5_,,$@) $@

rMatGraph_E_10_% : ../rMatGraph
	../rMatGraph  -a .5 -b .1 -m `expr 10 \* $(subst rMatGraph_E_10_,,$@)` $(subst rMatGraph_E_10_,,$@) $@

rMatGraph_E_12_% : ../rMatGraph
	../rMatGraph -a .55 -b .125 -m `expr 12 \* $(subst rMatGraph_E_12_,,$@)` $(subst rMatGraph_E_12_,,$@) $@

rMatGraph_E_16_% : ../rMatGraph
	../rMatGraph -a .57 -b .19 -m `expr 16 \* $(subst rMatGraph_E_16_,,$@)` $(subst rMatGraph_E_16_,,$@) $@

randLocalGraph_E_5_% : ../randLocalGraph
	../randLocalGraph -d 3 -m `expr 5 \* $(subst randLocalGraph_E_5_,,$@)` $(subst randLocalGraph_E_5_,,$@) $@

randLocalGraph_E_10_% : ../randLocalGraph
	../randLocalGraph -d 3 -m `expr 10 \* $(subst randLocalGraph_E_10_,,$@)` $(subst randLocalGraph_E_10_,,$@) $@

2Dgrid_WE_% : ../addWeights
	make -s $(subst WE,E,$@) 
	../addWeights $(subst WE,E,$@) $@

3Dgrid_WE_% : ../addWeights 
	make -s $(subst WE,E,$@) 
	../addWeights $(subst WE,E,$@) $@

rMatGraph_WE_% : ../addWeights
	make -s $(subst WE,E,$@) 
	../addWeights $(subst WE,E,$@) $@

randLocalGraph_WE_5_% : ../addWeights
	make -s $(subst WE,E,$@) 
	../addWeights $(subst WE,E,$@) $@

randLocalGraph_WE_10_% : ../addWeights
	make -s $(subst WE,E,$@) 
	../addWeights $(subst WE,E,$@) $@

2Dgrid_J_% : ../gridGraph
	../gridGraph -j -d 2 -o $(subst 2Dgrid_J_,,$@) $@

3Dgrid_J_% : ../gridGraph
	../gridGraph -j -d 3 -o $(subst 3Dgrid_J_,,$@) $@

3Dgrid_JR_% : ../gridGraph
	../gridGraph -j -d 3 $(subst 3Dgrid_JR_,,$@) $@

rMatGraph_J_5_% : ../rMatGraph
	../rMatGraph -j -a .5 -b .1 -m `expr 5 \* $(subst rMatGraph_J_5_,,$@)` $(subst rMatGraph_J_5_,,$@) $@

rMatGraph_J_12_% : ../rMatGraph
	../rMatGraph -j -a .55 -b .125 -o -m `expr 12 \* $(subst rMatGraph_J_12_,,$@)` $(subst rMatGraph_J_12_,,$@) $@

rMatGraph_JR_12_% : ../rMatGraph
	../rMatGraph -j -a .55 -b .125 -m `expr 12 \* $(subst rMatGraph_JR_12_,,$@)` $(subst rMatGraph_JR_12_,,$@) $@

rMatGraph_J_16_% : ../rMatGraph
	../rMatGraph -j -a .57 -b .19 -m `expr 16 \* $(subst rMatGraph_J_16_,,$@)` $(subst rMatGraph_J_16_,,$@) $@

randLocalGraph_J_5_% : ../randLocalGraph
	../randLocalGraph -j -d 3 -m `expr 5 \* $(subst randLocalGraph_J_5_,,$@)` $(subst randLocalGraph_J_5_,,$@) $@

randLocalGraph_J_10_% : ../randLocalGraph
	../randLocalGraph -j -o -d 3 -m `expr 10 \* $(subst randLocalGraph_J_10_,,$@)` $(subst randLocalGraph_J_10_,,$@) $@

randLocalGraph_JR_10_% : ../randLocalGraph
	../randLocalGraph -j -d 3 -m `expr 10 \* $(subst randLocalGraph_JR_10_,,$@)` $(subst randLocalGraph_JR_10_,,$@) $@

rMatGraph_J_500_% : ../rMatGraph
	../rMatGraph -j -a .5 -b .1 -m `expr 500 \* $(subst rMatGraph_J_500_,,$@)` $(subst rMatGraph_J_500_,,$@) $@

randLocalGraph_J_500_% : ../randLocalGraph
	../randLocalGraph -j -d 3 -m `expr 500 \* $(subst randLocalGraph_J_500_,,$@)` $(subst randLocalGraph_J_500_,,$@) $@

rMatGraph_E_500_% : ../rMatGraph
	../rMatGraph -a .5 -b .1 -m `expr 500 \* $(subst rMatGraph_E_500_,,$@)` $(subst rMatGraph_E_500_,,$@) $@

randLocalGraph_E_500_% : ../randLocalGraph
	../randLocalGraph -d 3 -m `expr 500 \* $(subst randLocalGraph_E_500_,,$@)` $(subst randLocalGraph_E_500_,,$@) $@

randDoubleVector_% : ../randDoubleVector
	../randDoubleVector $(subst randDoubleVector_,,$@) $@

starGraph_E_% : ../starGraph
	../starGraph $(subst starGraph_E_,,$@) $@

starGraph_WE_% : ../addWeights
	make -s $(subst WE,E,$@) 
	../addWeights $(subst WE,E,$@) $@

starGraph_J_% : ../starGraph
	../starGraph -j $(subst starGraph_J_,,$@) $@

starGraph_WJ_% : ../starGraph
	make -s $(subst WJ,J,$@)
	../adjGraphAddWeights $(subst WJ,J,$@) $@

combGraph_E_% : ../combGraph
	../combGraph -o $(subst combGraph_E_,,$@) $@

combGraph_J_K_% : ../combGraph
	../combGraph -o -j -k $(subst combGraph_J_K_,,$@) 25000000 $@

combGraph_WE_% : ../addWeights
	make -s $(subst WE,E,$@) 
	../addWeights $(subst WE,E,$@) $@

combGraph_J_% : ../combGraph
	../combGraph -o -j $(subst combGraph_J_,,$@) $@

expGraph_E_% : ../expGraph
	../expGraph $(subst expGraph_E_,,$@) $@

expGraph_WE_% : ../addWeights
	make -s $(subst WE,E,$@) 
	../addWeights $(subst WE,E,$@) $@

expGraph_J_% : ../expGraph
	../expGraph -j $(subst expGraph_J_,,$@) $@

randLocalGraph_WJ_5_% : ../adjGraphAddWeights
	make -s $(subst WJ,J,$@) 
	../adjGraphAddWeights $(subst WJ,J,$@) $@

rMatGraph_WJ_5_% : ../adjGraphAddWeights
	make -s $(subst WJ,J,$@) 
	../adjGraphAddWeights $(subst WJ,J,$@) $@

2Dgrid_WJ_% : ../adjGraphAddWeights
	make -s $(subst WJ,J,$@) 
	../adjGraphAddWeights $(subst WJ,J,$@) $@

3Dgrid_WJ_% : ../adjGraphAddWeights
	make -s $(subst WJ,J,$@) 
	../adjGraphAddWeights $(subst WJ,J,$@) $@

binTree_E_% : ../binTree
	../binTree $(subst binTree_E_,,$@) $@

binTree_WE_% : ../addWeights
	make -s $(subst WE,E,$@) 
	../addWeights $(subst WE,E,$@) $@

binTree_J_% : ../binTree
	../binTree -j $(subst binTree_J_,,$@) $@

binTree_WJ_% : ../adjGraphAddWeights
	make -s $(subst WJ,J,$@) 
	../adjGraphAddWeights $(subst WJ,J,$@) $@

randLocalGraph_F_% rMatGraph_F_%: ../adjGraphAddSourceSink
	make -s $(subst F,WJ,$@)
	../adjGraphAddSourceSink $(subst F,WJ,$@) $@

randLocalGraph_FD_% rMatGraph_FD_%: ../flowGraphToDimacs
	make -s $(subst FD,F,$@)
	../flowGraphToDimacs $(subst FD,F,$@) $@

ac_FD_%: ../maxFlowGens/bin/ac
	../maxFlowGens/bin/ac $(subst ac_FD_,,$@) 1000000 1 > $@

ak_FD_%: ../maxFlowGens/bin/ak
	../maxFlowGens/bin/ak $(subst ak_FD_,,$@) > $@

genrmf_long_FD_%: ../maxFlowGens/bin/genrmf-wrapper
	../maxFlowGens/bin/genrmf-wrapper long $(subst genrmf_long_FD_,,$@) > $@

genrmf_wide_FD_%: ../maxFlowGens/bin/genrmf-wrapper
	../maxFlowGens/bin/genrmf-wrapper wide $(subst genrmf_wide_FD_,,$@) > $@

washington_rlg_long_FD_%: ../maxFlowGens/bin/washington-wrapper
	../maxFlowGens/bin/washington-wrapper rlg-long $(subst washington_rlg_long_FD_,,$@) > $@

washington_rlg_wide_FD_%: ../maxFlowGens/bin/washington-wrapper
	../maxFlowGens/bin/washington-wrapper rlg-wide $(subst washington_rlg_wide_FD_,,$@) > $@

washington_line_moderate_FD_%: ../maxFlowGens/bin/washington-wrapper
	../maxFlowGens/bin/washington-wrapper line-moderate \
		$(subst washington_line_moderate_FD_,,$@) > $@

ac_F_% ak_F_% genrmf_long_F_% genrmf_wide_F_% \
	washington_rlg_long_F_% washington_rlg_wide_F_% \
	washington_line_moderate_F_%: ../dimacsToFlowGraph
	make -s $(subst F,FD,$@)
	../dimacsToFlowGraph $(subst F,FD,$@) $@
