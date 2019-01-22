/* maxflow generator in DIMACS format */
/*
   Implemented by 
   Tamas Badics, 1991, 
   Rutgers University, RUTCOR
   P.O.Box 5062
   New Brunswick, NJ, 08903
 
   e-mail: badics@rutcor.rutgers.edu
*/
/*
GENRMF -- Maxflow generator in DIMACS format.

Files: genio.c genrmf.c genmain.c  genio.h 
	   gen_maxflow_typedef.h  math_to_gcc.h
	   makefile

Compilation: Simply type make.

Usage: genrmf [-out out_file] [-seed seed]
              -a frame_size -b depth
              -c1 cap_range1 -c2 cap_range2

		Here without the -out option the generator will
		write to stdout.

		The generated network is as follows:
			It has b pieces of frames of size (a x a).
			(So alltogether the number of vertices is a*a*b)
			
			In each frame all the vertices are connected with 
			their neighbours. (forth and back)
			In addition the vertices of a frame are connected
			one to one with the vertices of next frame using 
			a random permutation of those vertices.

			The source is the lower left vertex of the first frame,
			the sink is the upper right vertex of the b'th frame. 
                         t
				+-------+
				|	   .|
				|	  .	|
             /  |    /  |
			+-------+/ -+ b
			|    |  |/.
		  a |   -v- |/
            |    |  |/
		    +-------+ 1
           s    a

			The capacities are randomly chosen integers
			from the range of (c1, c2) in the case of interconnecting
			edges, and c2 * a * a for the in-frame edges.
 
This generator was used by U. Derigs & W. Meier (1989)
in the article "Implementing Goldberg's Max-Flow-Algorithm
				A Computational Investigation"
ZOR - Methods & Models of OR (1989) 33:383-403

*/

#include <stdio.h>
#include "gen_maxflow_typedef.h"
#include "genio.h"
#include "math_to_gcc.h"

void make_edge(long from, long to, long c1, long c2);
void permute(void);
void connect(long offset, long cv, long x1, long y1);

network * N;
long * Parr;
long A, AA, C2AA, Ec;

/*==================================================================*/
network * gen_rmf(long a, long b, long c1, long c2)  
                              /* generates a network with 
								 a*a*b nodes and 6a*a*b-4ab-2a*a edges
								 random_frame network:
								 Derigs & Meier
								 Methods & Models of OR (1989) 
								 33:383-403 */
{
	long x, y, z, offset, cv;
	vertex * v;
	edge * e;
	
	A = a;
	AA = a*a;
	C2AA = c2*AA;
	Ec = 0;
	
	N = (network *)malloc(sizeof(network));
	N->vertnum = AA * b;
	N->edgenum = 5*AA*b-4*A*b-AA; 
	N->edges =(edge *)calloc(N->edgenum + 1, sizeof(edge));
	N->source = 1;
	N->sink   = N->vertnum;

	Parr = (long *)calloc(AA + 1, sizeof(long));
	
	for (x = 1; x <= AA; x++)
	  Parr[x] = x;
	
	for( z = 1; z <= b; z++){
		offset = AA * (z-1);
		if (z != b)
		  permute();
	
		for( x = 1; x <= A; x++){ 
			for( y = 1; y <= A; y++){
				cv = offset + (x - 1) * A + y;
				if (z != b)
				  make_edge(cv, offset + AA + Parr[cv - offset]  
							,c1, c2);   /* the intermediate edges */
				
				if (y < A)
				  connect(offset, cv, x, y + 1);
				if (y > 1)
				  connect(offset, cv, x, y - 1);
				if (x < A)
				  connect(offset, cv, x + 1, y);
				if (x > 1)
				  connect(offset, cv, x - 1, y);
			}
		}
	}
	return N;
}
/*==================================================================*/
void make_edge(long from, long to, long c1, long c2)
{
	Ec++;
	N->edges[Ec].from = from;
	N->edges[Ec].to = to;
	N->edges[Ec].cap = RANDOM(c1, c2);
}

/*==================================================================*/
void permute(void)
{
	long i, j, tmp;
	
	for (i = 1; i < AA; i++){
		j = RANDOM(i, AA);
		tmp = Parr[i];
		Parr[i] = Parr[j];
		Parr[j] = tmp;
	} 
}

/*==================================================================*/
void connect(long offset, long cv, long x1, long y1)
{
	long cv1;
	cv1 = offset + (x1 - 1) * A + y1;
	Ec++;
	N->edges[Ec].from = cv;
	N->edges[Ec].to = cv1;
	N->edges[Ec].cap = C2AA;
}


