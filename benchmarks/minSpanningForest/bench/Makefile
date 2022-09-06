include common/parallelDefs
BNCHMRK = MST

CHECKFILES = serialMST.o $(BNCHMRK)Check.o

COMMON = 

INCLUDE = -Icommon

%.o : %.C $(COMMON)
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

$(BNCHMRK)Check : $(CHECKFILES)
	$(CC) $(LFLAGS) -o $@ $(CHECKFILES)

clean :
	rm -f $(BNCHMRK)Check *.o *.pyc

