include common/parallelDefs

INCLUDE = -Icommon -Ipbbslib

wcCheck: wcCheck.C 
	$(CC) $(CFLAGS) $(LFLAGS) $(INCLUDE) -o wcCheck wcCheck.C

wcTime.o: wcTime.C
	$(CC) $(CFLAGS) $(INCLUDE) -o wcTime.o -c wcTime.C

clean :
	rm -f wcCheck *.o
