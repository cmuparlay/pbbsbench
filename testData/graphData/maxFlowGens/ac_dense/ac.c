/*
 ACYCLIC NETWORK GENERATOR FOR MAX-FLOW (by G. Waissi) 
                               (revised 11/25/90)
                               (revised 01/04/91)
			    (rewritten in C and modified by J. Setubal,
                           July 91. Program changed to generate instances
                           for DIMACS algorithm implementation challenge only.)

 usage: ac <num_nodes> <capacity> <seed>
 graph written to standard output
*/

#include <stdio.h>
#include <math.h>

int cap;
int p,q,i;
int num_arcs,num_nodes,capacity;
int tail,head;
int source;
int sink;
int input_seed;

main(argc,argv)
int argc;
char *argv[];
{

  num_nodes = atoi(argv[1]);
  capacity = atoi(argv[2]);
  input_seed = atoi(argv[3]);
  srandom(input_seed);
  UserValues();
}

Banner4()
{
  printf("c Fully Dense Acyclic Network\n");
  printf("c for Max-Flow\n");
  printf("c Arcs with random capacities [1:%d]\n", capacity);
  printf("c seed = %d\n", input_seed);
  printf("p max %d %d\n",num_nodes,num_arcs);
  printf("n %d s\n",source);
  printf("n %d t\n",sink);

}

AcyclicNet1()
{

  for (p = 1; p <= (num_nodes-1); p++)
    {
      tail = p;
      for (q = p+1; q <=  num_nodes; q++)
	{
	  head = q;
	  cap = RandomInteger(capacity);
	  printf("a %d %d %d\n",tail,head,cap);
	}
    }
}

UserValues()
{
  source = 1;
  sink = num_nodes;
  num_arcs = 0;
  for (i = 1; i <= (num_nodes-1); i++) num_arcs += i;
  Banner4();
  AcyclicNet1();
}

/* RandomInteger -- return a random integer from the range 1 .. high.
*/
int RandomInteger(high)
int high;
{
    return (random() % high) + 1;
}
