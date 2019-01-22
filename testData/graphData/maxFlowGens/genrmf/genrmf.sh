:	This is a shell archive.
:	Remove everything above this line and
:	run the following text with /bin/sh to create:
:	gen_maxflow_typedef.h
:	genio.c
:	genio.h
:	genmain.c
:	genrmf.c
:	index
:	makefile
:	math_to_gcc.h
:	readme
:	readme.shar
: This archive created: Sun Oct 13 17:26:13 1991
echo shar: extracting gen_maxflow_typedef.h
sed 's/^XX//' << 'SHAR_EOF' > gen_maxflow_typedef.h
XX/* gen_maxflow_typedef.h == Type definitions for a directed graph 
XX                            for generators */
XX/*
XX   Implemented by 
XX   Tamas Badics, 1991, 
XX   Rutgers University, RUTCOR
XX   P.O.Box 5062
XX   New Brunswick, NJ, 08903
XX 
XX   e-mail: badics@rutcor.rutgers.edu
XX*/
XX
XX#ifndef _GEN_MAXFLOW_TYPE_H
XX#define _GEN_MAXFLOW_TYPE_H
XX
XX#include <stdio.h>
XX
XX/*==================================================================*/
XXtypedef struct VERTEX{
XX	
XX	struct EDGE ** edgelist;  /* Pointer to the list of pointers to 
XX								 the adjacent edges. 
XX								 (No matter that to or from edges) */
XX
XX	struct EDGE ** current;   /* Pointer to the current edge */
XX
XX	int degree;        /* Number of adjacent edges (both direction) */
XX	int index;
XX}vertex;
XX
XX/*==================================================================*/
XXtypedef struct EDGE{
XX	int from;
XX	int to;
XX	int cap;        /* Capacity */
XX}edge;
XX
XX/*==================================================================*/
XXtypedef struct NETWORK{
XX
XX	struct NETWORK	* next, * prev;
XX
XX	int vertnum;
XX	int edgenum;
XX
XX	vertex	* verts; /* Vertex array[1..vertnum] */
XX	edge    * edges; /* Edge array[1..edgenum] */
XX
XX	int source; /* Pointer to the source */
XX	int sink;   /* Pointer to the sink */
XX}network;
XX
XX#endif 
XX
SHAR_EOF
if test 1250 -ne "`wc -c gen_maxflow_typedef.h`"
then
echo shar: error transmitting gen_maxflow_typedef.h '(should have been 1250 characters)'
fi
echo shar: extracting genio.c
sed 's/^XX//' << 'SHAR_EOF' > genio.c
XX/* I/O routines for DIMACS standard format generator */
XX/*
XX   Implemented by 
XX   Tamas Badics, 1991, 
XX   Rutgers University, RUTCOR
XX   P.O.Box 5062
XX   New Brunswick, NJ, 08903
XX 
XX   e-mail: badics@rutcor.rutgers.edu
XX*/
XX
XX#include <stdio.h>
XX#include "gen_maxflow_typedef.h"
XX#include "genio.h"
XX
XX/*===============================================================*/
XXvoid gen_free_net(network * n)
XX{
XX	free(n->edges);
XX	free(n);
XX}
XX
XX/*================================================================*/
XXvoid print_max_format (FILE * out, network * n
XX					   , char * comm[], int dim)
XX                           /* prints a network heading with 
XX							  dim lines of comments
XX							  (no \n needs at the ends )*/
XX	 
XX{
XX	int i, vnum, e_num;
XX	edge * e;
XX
XX	vnum = n->vertnum;
XX	e_num = n->edgenum;
XX	
XX	for( i = 0; i < dim; i++)
XX	  fprintf( out, "c %s\n", comm[i]);
XX	
XX	fprintf( out, "p max %7d %10d\n", vnum, e_num);
XX	fprintf( out, "n %7d s\n", n->source);
XX	fprintf( out, "n %7d t\n", n->sink);
XX
XX	for (i = 1; i <= e_num; i++){
XX		e = &n->edges[i];
XX		fprintf(out, "a %7d %7d %10d\n"
XX			   , e->from, e->to, (int)e->cap); 
XX	}
XX}
XX
SHAR_EOF
if test 1107 -ne "`wc -c genio.c`"
then
echo shar: error transmitting genio.c '(should have been 1107 characters)'
fi
echo shar: extracting genio.h
sed 's/^XX//' << 'SHAR_EOF' > genio.h
XX#ifndef _GENIO
XX#define _GENIO
XX
XXvoid gen_free_net(network * n);
XX
XXvoid print_max_format (FILE * outfile, network * n
XX					   , char * comm[], int dim);
XX                           /* prints a network heading with 
XX							  dim lines of comments
XX							  (no \n needs at the ends )*/
XX	 
XXnetwork * gen_rmf(int a, int b, int c1, int c2);  
XX                              /* generates a network with 
XX								 a*a*b nodes and 6a*a*b-4ab-2a*a edges
XX								 random_frame network:
XX								 Derigs & Meier
XX								 Methods & Models of OR (1989) 
XX								 33:383-403 */
XX#endif
SHAR_EOF
if test 565 -ne "`wc -c genio.h`"
then
echo shar: error transmitting genio.h '(should have been 565 characters)'
fi
echo shar: extracting genmain.c
sed 's/^XX//' << 'SHAR_EOF' > genmain.c
XX/* genrmf maxflow input generator */
XX
XX#include <stdio.h>
XX#include "gen_maxflow_typedef.h"
XX#include "genio.h"
XX
XX
XXvoid print_usage(void);
XX
XXvoid main(int argc, char * argv[])
XX{
XX	network * n;
XX	int i, feas, quiet;
XX	FILE * output;
XX	int a, b, c1, c2;
XX	char comm[10][80];
XX	char * com1[10];
XX	int seed;
XX
XX	output = stdout;
XX	a = -1;
XX	b = -1;
XX	c1 = -1;
XX	c2 = -1;
XX	seed = -1;
XX	
XX	for (i = 1; i < argc; i++){
XX		switch (look_up(argv[i])){
XX		  case 0: 
XX			output = fopen(argv[++i],"w");
XX			if (output == NULL) {
XX				fprintf(stderr
XX                    ,"genrmf: Output file %s can't be opened\n",argv[i]);
XX				exit(-1);
XX			}	
XX			break;
XX		  case 1: 
XX			a = atoi(argv[++i]);
XX			break;
XX		  case 2: 
XX			b = atoi(argv[++i]);
XX			break;
XX		  case 3: 
XX			c1 = atoi(argv[++i]);
XX			break;
XX		  case 4: 
XX			c2 = atoi(argv[++i]);
XX			break;
XX		  case 5: 
XX			seed = atoi(argv[++i]);
XX			break;
XX		  default:
XX			break;
XX		} 
XX	}
XX
XX	if (a == -1 || b == -1 || c1 == -1 || c2 == -1)
XX	  print_usage();
XX	
XX	if (seed = -1)
XX	  seed = (int) time(0);
XX	
XX	srand48(seed);
XX	
XX	n = gen_rmf(a, b, c1, c2);
XX
XX	sprintf(comm[0], "This file was generated by genrmf.");
XX	sprintf(comm[1], "The parameters are: a: %d b: %d c1: %d c2: %d"
XX			, a, b, c1, c2);
XX
XX	com1[0] = comm[0];
XX	com1[1] = comm[1];
XX	
XX	print_max_format(output, n, com1, 2);
XX
XX	gen_free_net(n);
XX
XX	if (output != stdout)
XX	  fclose(output);
XX}
XX/*=================================================================*/
XX#define OPS_NUM 6
XX
XXint look_up(char * s)
XX{
XX	char * ops[OPS_NUM] 
XX	  = { "-out", "-a", "-b", "-c1", "-c2", "-seed"};
XX	int i;
XX	
XX	for (i = 0; i < OPS_NUM; i++){
XX		if (strcmp(ops[i], s) == 0)
XX		  return i;
XX	} 
XX	return -1;
XX} 
XX
XXvoid print_usage(void)
XX{
XX	printf("Usage: genrmf [-out out_file]\n");
XX	printf("              -a frame_size -b depth\n");
XX	printf("              -c1 cap_range1 -c2 cap_range2\n");
XX	exit(0);
XX}
XX
SHAR_EOF
if test 1819 -ne "`wc -c genmain.c`"
then
echo shar: error transmitting genmain.c '(should have been 1819 characters)'
fi
echo shar: extracting genrmf.c
sed 's/^XX//' << 'SHAR_EOF' > genrmf.c
XX/* maxflow generator in DIMACS format */
XX/*
XX   Implemented by 
XX   Tamas Badics, 1991, 
XX   Rutgers University, RUTCOR
XX   P.O.Box 5062
XX   New Brunswick, NJ, 08903
XX 
XX   e-mail: badics@rutcor.rutgers.edu
XX*/
XX/*
XXGENRMF -- Maxflow generator in DIMACS format.
XX
XXFiles: genio.c genrmf.c genmain.c  genio.h 
XX	   gen_maxflow_typedef.h  math_to_gcc.h
XX	   makefile
XX
XXCompilation: Simply type make.
XX
XXUsage: genrmf [-out out_file]
XX              -a frame_size -b depth
XX              -c1 cap_range1 -c2 cap_range2
XX
XX		Here without the -out option the generator will
XX		write to stdout.
XX
XX		The generated network is as follows:
XX			It has b pieces of frames of size (a x a).
XX			(So alltogether the number of vertices is a*a*b)
XX			
XX			In each frame all the vertices are connected with 
XX			their neighbours. (forth and back)
XX			In addition the vertices of a frame are connected
XX			one to one with the vertices of next frame using 
XX			a random permutation of those vertices.
XX
XX			The source is the lower left vertex of the first frame,
XX			the sink is the upper right vertex of the b'th frame. 
XX                         t
XX				+-------+
XX				|	   .|
XX				|	  .	|
XX             /  |    /  |
XX			+-------+/ -+ b
XX			|    |  |/.
XX		  a |   -v- |/
XX            |    |  |/
XX		    +-------+ 1
XX           s    a
XX
XX			The capacities are randomly chosen integers
XX			from the range of (c1, c2) in the case of interconnecting
XX			edges, and c2 * a * a for the in-frame edges.
XX 
XXThis generator was used by U. Derigs & W. Meier (1989)
XXin the article "Implementing Goldberg's Max-Flow-Algorithm
XX				A Computational Investigation"
XXZOR - Methods & Models of OR (1989) 33:383-403
XX
XX*/
XX
XX#include <stdio.h>
XX#include "gen_maxflow_typedef.h"
XX#include "genio.h"
XX#include "math_to_gcc.h"
XX
XXvoid make_edge(int from, int to, int c1, int c2);
XXvoid permute(void);
XXvoid connect(int offset, int cv, int x1, int y1);
XX
XXnetwork * N;
XXint * Parr;
XXint A, AA, C2AA, Ec;
XX
XX/*==================================================================*/
XXnetwork * gen_rmf(int a, int b, int c1, int c2)  
XX                              /* generates a network with 
XX								 a*a*b nodes and 6a*a*b-4ab-2a*a edges
XX								 random_frame network:
XX								 Derigs & Meier
XX								 Methods & Models of OR (1989) 
XX								 33:383-403 */
XX{
XX	int x, y, z, offset, cv;
XX	vertex * v;
XX	edge * e;
XX	
XX	A = a;
XX	AA = a*a;
XX	C2AA = c2*AA;
XX	Ec = 0;
XX	
XX	N = (network *)malloc(sizeof(network));
XX	N->vertnum = AA * b;
XX	N->edgenum = 5*AA*b-4*A*b-AA; 
XX	N->edges =(edge *)calloc(N->edgenum + 1, sizeof(edge));
XX	N->source = 1;
XX	N->sink   = N->vertnum;
XX
XX	Parr = (int *)calloc(AA + 1, sizeof(int));
XX	
XX	for (x = 1; x <= AA; x++)
XX	  Parr[x] = x;
XX	
XX	for( z = 1; z <= b; z++){
XX		offset = AA * (z-1);
XX		if (z != b)
XX		  permute();
XX	
XX		for( x = 1; x <= A; x++){ 
XX			for( y = 1; y <= A; y++){
XX				cv = offset + (x - 1) * A + y;
XX				if (z != b)
XX				  make_edge(cv, offset + AA + Parr[cv - offset]  
XX							,c1, c2);   /* the intermediate edges */
XX				
XX				if (y < A)
XX				  connect(offset, cv, x, y + 1);
XX				if (y > 1)
XX				  connect(offset, cv, x, y - 1);
XX				if (x < A)
XX				  connect(offset, cv, x + 1, y);
XX				if (x > 1)
XX				  connect(offset, cv, x - 1, y);
XX			}
XX		}
XX	}
XX	return N;
XX}
XX/*==================================================================*/
XXvoid make_edge(int from, int to, int c1, int c2)
XX{
XX	Ec++;
XX	N->edges[Ec].from = from;
XX	N->edges[Ec].to = to;
XX	N->edges[Ec].cap = RANDOM(c1, c2);
XX}
XX
XX/*==================================================================*/
XXvoid permute(void)
XX{
XX	int i, j, tmp;
XX	
XX	for (i = 1; i < AA; i++){
XX		j = RANDOM(i, AA);
XX		tmp = Parr[i];
XX		Parr[i] = Parr[j];
XX		Parr[j] = tmp;
XX	} 
XX}
XX
XX/*==================================================================*/
XXvoid connect(int offset, int cv, int x1, int y1)
XX{
XX	int cv1;
XX	cv1 = offset + (x1 - 1) * A + y1;
XX	Ec++;
XX	N->edges[Ec].from = cv;
XX	N->edges[Ec].to = cv1;
XX	N->edges[Ec].cap = C2AA;
XX}
XX
XX
SHAR_EOF
if test 3843 -ne "`wc -c genrmf.c`"
then
echo shar: error transmitting genrmf.c '(should have been 3843 characters)'
fi
echo shar: extracting index
sed 's/^XX//' << 'SHAR_EOF' > index
XXtotal 31
XX-rw-rw-r--  1 badics       1250 May 29 19:55 gen_maxflow_typedef.h
XX-rw-rw-r--  1 badics       1107 May 29 19:55 genio.c
XX-rw-rw-r--  1 badics        565 May 29 19:55 genio.h
XX-rw-rw-r--  1 badics       1676 May 29 19:55 genmain.c
XX-rw-rw-r--  1 badics       3843 May 29 19:55 genrmf.c
XX-rw-rw-r--  1 mcgeoch     13153 Sep 13 10:25 genrmf.sh
XX-r--rw-r--  1 badics          0 Oct 13 13:01 index
XX-rw-rw-r--  1 badics        189 May 29 19:55 makefile
XX-rw-rw-r--  1 badics        318 May 29 19:55 math_to_gcc.h
XX-rw-rw-r--  1 mcgeoch      1371 Jul  1 19:13 readme
XX-rw-rw-r--  1 mcgeoch      1418 Jul  1 19:04 readme.BAK
XX-rw-rw-r--  1 mcgeoch       198 Sep 13 10:29 readme.shar
SHAR_EOF
if test 675 -ne "`wc -c index`"
then
echo shar: error transmitting index '(should have been 675 characters)'
fi
echo shar: extracting makefile
sed 's/^XX//' << 'SHAR_EOF' > makefile
XX
XXCFILES1= genio.c genrmf.c genmain.c
XX
XXOBJS= $(CFILES1:%.c=%.o)
XX
XXCC= gcc
XX
XXOUTFILE = genrmf
XX
XXCFLAGS=-O 
XX
XXLIBS = -lm 
XX
XXmain: ${OBJS} 
XX	${CC} ${OBJS} ${CFLAGS} ${LDFLAGS} -o ${OUTFILE} ${LIBS}
SHAR_EOF
if test 189 -ne "`wc -c makefile`"
then
echo shar: error transmitting makefile '(should have been 189 characters)'
fi
echo shar: extracting math_to_gcc.h
sed 's/^XX//' << 'SHAR_EOF' > math_to_gcc.h
XX/* Macros for converting from Turboc math.h to gcc math.h */
XX
XX#ifndef MATH_TO
XX#define MATH_TO
XX
XX#include <math.h>
XX
XXdouble drand48(void);
XXvoid srand48(long int);
XX
XX#define random(A) (int)(drand48()*(double)(A))
XX#define RANDOM(A, B) (int)(random(B-A + 1) + (A))
XX
XX#define sgn(A)  ((A) > 0) ? 1 : ((A)== 0) ? 0 : -1
XX
XX#endif
SHAR_EOF
if test 318 -ne "`wc -c math_to_gcc.h`"
then
echo shar: error transmitting math_to_gcc.h '(should have been 318 characters)'
fi
echo shar: extracting readme
sed 's/^XX//' << 'SHAR_EOF' > readme
XXGENRMF -- Maxflow generator in DIMACS format.
XX
XXThis generator produces the RMFGEN networks developed by 
XX Goldfarb and Grigoriadis (see ``A computational comparison of the Dinic
XX  and Network Simplex methods for maximum flow,'' Annals of Operations 
XX  Research 13 (1988), pp 83--123. 
XX
XXContributed by Tamas Badics. 
XX
XXFiles: genio.c genrmf.c genmain.c  genio.h 
XX	   gen_maxflow_typedef.h  math_to_gcc.h
XX	   makefile
XX
XXCompilation: Get all the files. Type make.
XX
XXUsage: genrmf [-out out_file]
XX              -a frame_size -b depth
XX              -c1 cap_range1 -c2 cap_range2
XX
XX	       Here without the -out option the generator will
XX		write to stdout.
XX
XX		The generated network is as follows:
XX			It has b pieces of frames of size (a x a).
XX			(So alltogether the number of vertices is a*a*b)
XX			
XX			In each frame all the vertices are connected with 
XX			their neighbours. (forth and back)
XX			In addition the vertices of a frame are connected
XX			one to one with the vertices of next frame using 
XX			a random permutation of those vertices.
XX
XX			The source is the lower left vertex of the first frame,
XX			the sink is the upper right vertex of the b'th frame. 
XX               
XX			The capacities are randomly chosen integers
XX			from the range of (c1, c2) in the case of 
XX                        interconnecting edges, and c2 * a * a for 
XX                        the in-frame edges.
XX
XX
XX
SHAR_EOF
if test 1371 -ne "`wc -c readme`"
then
echo shar: error transmitting readme '(should have been 1371 characters)'
fi
echo shar: extracting readme.shar
sed 's/^XX//' << 'SHAR_EOF' > readme.shar
XX
XXThe file genrmf.sh is a shar file containing all the files
XXin this directory.  You can get the file, then type
XX
XXmkdir genrfm
XX/bin/sh genrmf.sh
XX
XXTo create a subdirectory containing all the files. 
XX
SHAR_EOF
if test 198 -ne "`wc -c readme.shar`"
then
echo shar: error transmitting readme.shar '(should have been 198 characters)'
fi
:	End of shell archive
exit 0
