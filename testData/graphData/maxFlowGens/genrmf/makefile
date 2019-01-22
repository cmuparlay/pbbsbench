
CFILES1= genio.c genrmf.c genmain.c

OBJS= $(CFILES1:%.c=%.o)

OUTFILE = genrmf

CFLAGS=-O -g

LIBS = -lm 

main: ${OBJS} 
	${CC} ${OBJS} ${CFLAGS} ${LDFLAGS} -o ${OUTFILE} ${LIBS}
clean:
	rm ${OUTFILE} ${OBJS}