#ifndef _GENIO
#define _GENIO

void gen_free_net(network * n);

void print_max_format (FILE * outfile, network * n
					   , char * comm[], int dim);
                           /* prints a network heading with 
							  dim lines of comments
							  (no \n needs at the ends )*/
	 
network * gen_rmf(long a, long b, long c1, long c2);  
                              /* generates a network with 
								 a*a*b nodes and 6a*a*b-4ab-2a*a edges
								 random_frame network:
								 Derigs & Meier
								 Methods & Models of OR (1989) 
								 33:383-403 */
#endif
