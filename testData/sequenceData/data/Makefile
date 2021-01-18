GENERATORS = ../randomSeq ../equalSeq ../almostEqualSeq ../almostSortedSeq ../exptSeq ../trigramSeq ../addDataSeq ../trigramString

STRINGFILES = wikipedia250M.txt wikisamp.xml chr22.dna etext99 
STRINGFILES_LONG = wikisamp.xml chr22.dna etext99 HG18 howto jdk13c proteins rctail96 rfc sprot34 w3c2
STRINGBZIP = $(addsuffix .bz2,$(STRINGFILES))

CLASSIFY = covtype.data kddcup.data 
CLASSIFYFILES = $(addsuffix .train,$(CLASSIFY)) $(addsuffix .test,$(CLASSIFY)) $(addsuffix .labels,$(CLASSIFY))
CLASSIFYBZIP = $(addsuffix .bz2,$(CLASSIFYFILES))

$(GENERATORS) :
	cd ..; make -s $(subst ../,,$@)

$(STRINGFILES) : $(STRINGBZIP)
	bunzip2 -k -f $@.bz2
	touch $@

$(STRINGBZIP) :
	cp ../../data/$@ .

%.data : %.data.train %.data.test %.data.labels $(CLASSIFYBZIP)
	echo

$(CLASSIFYFILES) : $(CLASSIFYBZIP)
	bunzip2 -k -f $@.bz2
	touch $@

$(CLASSIFYBZIP) :
	cp ../../data/$@ .

randomSeq_%_pair_int :  ../addDataSeq
	make -s $(subst _pair_int,,$@)
	../addDataSeq -t int $(subst _pair_int,,$@) $@

randomSeq_%_pair_double :  ../addDataSeq
	make -s $(subst _pair_double,,$@)
	../addDataSeq -t double $(subst _pair_double,,$@) $@

randomSeq_10_% : ../randomSeq
	../randomSeq -t $(subst randomSeq_10_,,$@) 10 $@

randomSeq_100_% : ../randomSeq
	../randomSeq -t $(subst randomSeq_100_,,$@) 100 $@

randomSeq_1000_% : ../randomSeq
	../randomSeq -t $(subst randomSeq_1000_,,$@) 1000 $@

randomSeq_10000_% : ../randomSeq
	../randomSeq -t $(subst randomSeq_10000_,,$@) 10000 $@

randomSeq_100000_% : ../randomSeq
	../randomSeq -t $(subst randomSeq_100000_,,$@) 100000 $@

randomSeq_1000000_% : ../randomSeq
	../randomSeq -t $(subst randomSeq_1000000_,,$@) 1000000 $@

randomSeq_10M_% : ../randomSeq
	../randomSeq -t $(subst randomSeq_10M_,,$@) 10000000 $@

randomSeq_100M_% : ../randomSeq
	../randomSeq -t $(subst randomSeq_100M_,,$@) 100000000 $@

randomSeq_% : ../randomSeq
	../randomSeq -t int $(subst randomSeq_,,$@) $@

equalSeq_% : ../equalSeq
	../equalSeq -t int -r $(subst equalSeq_,,$@) $(subst equalSeq_,,$@) $@

almostEqualSeq_% : ../almostEqualSeq
	../almostEqualSeq -t int -r $(subst almostEqualSeq_,,$@) $(subst almostEqualSeq_,,$@) $@

# In this one the integers are limited to a range of 100K
randomSeq_10M_100K_int : ../randomSeq
	../randomSeq -t int -r 100000 10000000 $@

randomSeq_100M_100K_int : ../randomSeq
	../randomSeq -t int -r 100000 100000000 $@

randomSeq_10M_256_int : ../randomSeq
	../randomSeq -t int -r 256 10000000 $@

randomSeq_100M_256_int : ../randomSeq
	../randomSeq -t int -r 256 100000000 $@

exptSeq_10M_% : ../exptSeq
	../exptSeq -t $(subst exptSeq_10M_,,$@) 10000000 $@

exptSeq_100M_% : ../exptSeq
	../exptSeq -t $(subst exptSeq_100M_,,$@) 100000000 $@

almostSortedSeq_10M_% : ../almostSortedSeq
	../almostSortedSeq -t $(subst almostSortedSeq_10M_,,$@) 10000000 $@

almostSortedSeq_100M_% : ../almostSortedSeq
	../almostSortedSeq -t $(subst almostSortedSeq_100M_,,$@) 100000000 $@

trigramSeq_10M : ../trigramSeq
	../trigramSeq 10000000 $@

trigramSeq_100M : ../trigramSeq
	../trigramSeq 100000000 $@

trigramSeq_% : ../trigramSeq
	../trigramSeq $(subst trigramSeq_,,$@) $@

trigramString_100M : ../trigramString
	../trigramString 100000000 $@

trigramString_1G : ../trigramString
	../trigramString 1000000000 $@

trigramString_% : ../trigramString
	../trigramString $(subst trigramString_,,$@) $@

clean :
	rm -f *0* $(STRINGFILES) $(STRINGBZIP) $(CLASSIFYFILES) $(CLASSIFYBZIP)
