include common/parallelDefs

COMMON = common/IO.h common/parse_command_line.h common/geometry.h common/geometryIO.h parlay/parallel.h parlay/primitives.h
GENERATORS = randPoints
TO_UPDATE = triangles addRays toNodes

.PHONY: all clean
all: $(GENERATORS)

%.o : %.C $(COMMON)
	$(CC) $(CFLAGS) -c $< -o $@

randPoints : randPoints.C geometryData.h $(COMMON)
	$(CC) $(CFLAGS) -o $@ $@.C $(LFLAGS)

toNodes : toNodes.C $(COMMON)
	$(CC) $(CFLAGS) -o $@ $@.C $(LFLAGS)

triangles : triangles.C geometryData.h $(COMMON)
	$(CC) $(CFLAGS) -o $@ $@.C $(LFLAGS)

addRays : addRays.C $(COMMON)
	$(CC) $(CFLAGS) -o $@ $@.C $(LFLAGS)

clean :
	rm -f *.o $(GENERATORS)
	make clean -s -C data

cleansrc :
	make -s clean
	rm -f $(COMMON)
