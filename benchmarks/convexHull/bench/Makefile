include common/parallelDefs
BNCHMRK = hull

CHECKFILES = $(BNCHMRK)Check.o

COMMON = 

INCLUDE = -Icommon

%.o : %.C $(COMMON)
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

$(BNCHMRK)Check : $(CHECKFILES)
	$(CC) -o $@ $(CHECKFILES) $(LFLAGS)

clean :
	rm -f $(BNCHMRK)Check *.o *.pyc
