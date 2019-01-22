
/*

  June 29, 1999 -- Andrew Goldberg.
    Removed MAX_N parameter - not needed with dynamic allocation.

    Changed MAX_DEGREE to a variable deg, used for dynamic memory
    allocation.

 */

/*
Ver 2.2
Bug fixes - basic line graph

Ver 2.1

Added myAlloc to check for allocation failure

Ver 2.0
 
Jun 28, 1999 -- Richard Anderson
  Converted from static to dynamic memory allocation of graph structures.

  The graph data structure is an array for the vertices and an 
  array for the adjacency lists.

  A number of unused routines were removed from the file.  (They would
  have needed to be rewritten to support the changes).  

  I did not attempt to document, or otherwise clean up this ancient code.
  */


  /* makegraph.c */

/*
  The seed input parameter added by Cherkassy and Goldberg for reproducibility
 */

/* graph.h */

#include <sys/time.h>
#include <stdio.h>
/*#include <time.h>
*/
#include <ctype.h>
#include <strings.h>

/* ghead.h */



#define FAILURE    0
#define SUCCESS    1
#define FALSE      0
#define TRUE       1



#define MAX_CAP   100000000

/* Dimacs problem types */
#define UNDEFINED        0
#define MINCOSTFLOW      1
#define MAXFLOW          2
#define ASSIGNMENT       3

typedef struct enode {
  struct enode *next;
  struct enode *mate;
  int c;
  int f;
  int h;
  int t;
  int flag;
} Edge;



typedef struct {
  Edge **A;
  int *V;
  int size;
  int max_v;
  int max_n;           /* size of arrays allocated for vertices  */
                       /* and adjacency lists */
} Graph;

typedef struct {
  int head, tail, size;
  int *data;
} Queue;





#define MAX_RANGE 1000000

#define VERY_BIG 10000000

#define LOG_RANGE 20             /* Number of discrete values available */

int Range[] = {1000000, 500000, 250000, 125000, 62500, 31250,
		 15625, 7812, 3906, 1953, 976, 488, 244, 122,
		 61, 31, 15, 7, 4, 2};


void *myAlloc();

main(argc, argv)
int argc;
char *argv[];
{
  Graph *G, *Mesh(), *RLevel(), *R2Level(), *Match(), *SquareMesh(), 
        *BasicLine(), *ExponentialLine(), *DExponentialLine(),
*DinicBadCase(),
        *GoldBadCase(), *Cheryian();

  FILE *f;
  int dim1, dim2, range, fct, s, t;
  int init_seed;
  void GraphOutput();

  if (argc != 6)
    {
      fprintf( stderr, "Usage: %s fct dim1 dim2 range seed\n", argv[0]);
      exit (2);
    }

  fct = atoi(argv[1]);
  dim1 = atoi(argv[2]);
  dim2 = atoi(argv[3]);
  range = atoi(argv[4]);
  init_seed = atoi(argv[5]);
  if ( init_seed < 0 ) init_seed = -init_seed;
  init_seed = 2 * init_seed + 1;


  InitRandom(init_seed);

  f = stdout;

  switch(fct){
  case 1:
    fprintf(f, "c Mesh Graph\n");
    fprintf(f, "c %d Rows, %d columns, capacities in range [0, %d]\n",
	    dim1, dim2, range);
    G = Mesh(dim1, dim2, range);
    s = 0;
    t = G->size - 1;
    break;
  case 2:
    fprintf(f, "c Random Leveled Graph\n");
    fprintf(f, "c %d Rows, %d columns, capacities in range [0, %d]\n",
	    dim1, dim2, range);
    G = RLevel(dim1, dim2, range);
    s = 0;
    t = G->size - 1;
    break;
  case 3:
    fprintf(f, "c Random 2 Leveled Graph\n");
    fprintf(f, "c %d Rows, %d columns, capacities in range [0, %d]\n",
	    dim1, dim2, range);
    G = RLevel(dim1, dim2, range);
    s = 0;
    t = G->size - 1;
    break;
  case 4:
    fprintf(f, "c Matching Graph\n");
    fprintf(f, "c %d vertices, %d degree\n",
	    dim1, dim2);
    G = Match(dim1, dim2);
    s = 0;
    t = G->size - 1;
    break;

  case 5:
    fprintf(f, "c Square Mesh\n");
    fprintf(f, "c %d x %d vertices, %d degree, range [0,%d]\n", 
	    dim1, dim1, dim2, range);
    G = SquareMesh(dim1, dim2, range);
    s = 0;
    t = G->size - 1;
    break;

  case 6:
    fprintf(f, "c Basic Line Mesh\n");
    fprintf(f, "c %d x %d vertices, degree d\n", 
	    dim1, dim2, range);
    G = BasicLine(dim1, dim2, range);
    s = 0;
    t = G->size - 1;
    break;

  case 7:
    fprintf(f, "c Exponential Line\n");
    fprintf(f, "c %d x %d vertices, degree %d\n", 
	    dim1, dim2, range);
    G = ExponentialLine(dim1, dim2, range);
    s = 0;
    t = G->size - 1;
    break;

  case 8:
    fprintf(f, "c Double Exponential Line\n");
    fprintf(f, "c %d x %d vertices, degree %d\n", 
	    dim1, dim2, range);
    G = DExponentialLine(dim1, dim2, range);
    s = 0;
    t = G->size - 1;
    break;

  case 9:
    fprintf(f, "c Line Graph - Bad case for Dinics\n");
    fprintf(f, "c %d vertices\n", dim1);
    G = DinicBadCase(dim1);
    s = 0;
    t = G->size - 1;
    break;

  case 10:
    fprintf(f, "c  Bad case for Goldberg\n");
    fprintf(f, "c %d vertices\n", dim1);
    G = GoldBadCase(dim1);
    s = 0;
    t = G->size - 1;
    break;

  case 11:
    fprintf(f, "c  Cheryian Graph\n");
    fprintf(f, "c n = %d, m = %d, c = %d, total vertices %d \n", 
	    dim1, dim2, range, 4*dim2*range + 6 + dim1);
    G = Cheryian(dim1, dim2, range);    
    s = 0;
    t = G->size - 1;
    break;

  default:
    Barf("Undefined class");
    break;

  }

  GraphOutput(G, f, s, t);
}

Graph *Mesh(d1, d2, r)
int d1, d2, r;
{
  Graph *G;
  int i, j, source, sink;

  if (d1 < 2 || d2 < 2)
    Barf("Degenerate graph");

  G = (Graph *) myAlloc(sizeof(Graph));

  InitGraph(G, d1*d2 + 2);  

  for (i = 0; i <= d1*d2 + 1; i++)
    AddVertex(i, G);

  source = 0;
  sink = d1*d2 + 1;

  for (i = 1; i <= d1; i++){
    AddEdge(source, source + i, 3*r, G);
    AddEdge(sink - i, sink, 3*r, G);
  }

  for (i = 0; i < d2 - 1; i++){
    AddEdge(i*d1 + 1, (i+1)*d1 + d1, RandomInteger(1, r), G);
    AddEdge(i*d1 + 1, (i+1)*d1 + 1, RandomInteger(1, r), G);
    AddEdge(i*d1 + 1, (i+1)*d1 + 2, RandomInteger(1, r), G);
    for (j = 2; j <= d1 - 1; j++){
      AddEdge(i*d1 + j, (i+1)*d1 + j - 1, RandomInteger(1, r), G);
      AddEdge(i*d1 + j, (i+1)*d1 + j, RandomInteger(1, r), G);
      AddEdge(i*d1 + j, (i+1)*d1 + j + 1, RandomInteger(1, r), G);
    }
    AddEdge(i*d1 + d1, (i+1)*d1 + d1 - 1, RandomInteger(1, r), G);
    AddEdge(i*d1 + d1, (i+1)*d1 + d1, RandomInteger(1, r), G);
    AddEdge(i*d1 + d1, (i+1)*d1 + 1, RandomInteger(1, r), G);
  }

  return G;
}


Graph *RLevel(d1, d2, r)
int d1, d2, r;
{
  Graph *G;
  int i, j, source, sink, x[3];

  if (d1 < 2 || d2 < 2)
    Barf("Degenerate graph");

  G = (Graph *) myAlloc(sizeof(Graph));

  InitGraph(G, d1*d2 + 2);  

  for (i = 0; i <= d1*d2 + 1; i++)
    AddVertex(i, G);

  source = 0;
  sink = d1*d2 + 1;

  for (i = 1; i <= d1; i++){
    AddEdge(source, source + i, 3*r, G);
    AddEdge(sink - i, sink, 3*r, G);
  }

  for (i = 0; i < d2 - 1; i++){
    for (j = 1; j <= d1; j++){    
      RandomSubset(1, d1, 3, x);
      AddEdge(i*d1 + j, (i+1)*d1 + x[0], RandomInteger(1, r), G);
      AddEdge(i*d1 + j, (i+1)*d1 + x[1], RandomInteger(1, r), G);
      AddEdge(i*d1 + j, (i+1)*d1 + x[2], RandomInteger(1, r), G);
    }
  }

  return G;
}


Graph *R2Level(d1, d2, r)
int d1, d2, r;
{
  Graph *G;
  int i, j, source, sink, x[3];

  if (d1 < 2 || d2 < 2)
    Barf("Degenerate graph");

  G = (Graph *) myAlloc(sizeof(Graph));

  InitGraph(G, d1*d2+2);  

  for (i = 0; i <= d1*d2 + 1; i++)
    AddVertex(i, G);

  source = 0;
  sink = d1*d2 + 1;

  for (i = 1; i <= d1; i++){
    AddEdge(source, source + i, 3*r, G);
    AddEdge(sink - i, sink, 3*r, G);
  }

  for (i = 0; i < d2 - 2; i++){
    for (j = 1; j <= d1; j++){    
      RandomSubset(1, 2*d1, 3, x);
      AddEdge(i*d1 + j, (i+1)*d1 + x[0], RandomInteger(1, r), G);
      AddEdge(i*d1 + j, (i+1)*d1 + x[1], RandomInteger(1, r), G);
      AddEdge(i*d1 + j, (i+1)*d1 + x[2], RandomInteger(1, r), G);
    }
  }
  for (j = 1; j <= d1; j++){    
    RandomSubset(1, d1, 3, x);
      AddEdge((d2-2)*d1 + j, (d1-1)*d1 + x[0], RandomInteger(1, r), G);
      AddEdge((d2-2)*d1 + j, (d1-1)*d1 + x[1], RandomInteger(1, r), G);
      AddEdge((d2-2)*d1 + j, (d1-1)*d1 + x[2], RandomInteger(1, r), G);
    }

  return G;
}




Graph *Match(n, d)
int n, d;
{
  Graph *G;
  int i, j, source, sink;
  int *x;

  if (n < 2 || d > n)
    Barf("Degenerate graph");

  G = (Graph *) myAlloc(sizeof(Graph));

  InitGraph(G, 2*n + 2);  

  for (i = 0; i <= 2*n + 1; i++)
    AddVertex(i, G);

  source = 0;
  sink = 2*n + 1;

  for (i = 1; i <= n; i++){
    AddEdge(source, source + i, 1, G);
    AddEdge(sink - i, sink, 1, G);
  }


  x = (int *) myAlloc ((2*n + 2) * sizeof(int));

  for (j = 1; j <= n; j++){    
    RandomSubset(1, n, d, x);
    for (i = 0; i < d; i++)
      AddEdge(j, n + x[i], 1, G);
  }
  return G;
}


Graph *SquareMesh(d, deg, r)
int d, deg, r;
{
  Graph *G;
  int i, j, k, source, sink;

  if (d < deg)
    Barf("Degenerate graph");

  G = (Graph *) myAlloc(sizeof(Graph));

  InitGraph(G, d*d + 2);  

  for (i = 0; i <= d*d + 1; i++)
    AddVertex(i, G);

  source = 0;
  sink = d*d + 1;

  for (i = 1; i <= d; i++){
    AddEdge(source, source + i, 3*r, G);
    AddEdge(sink - i, sink, 3*r, G);
  }

  for (i = 0; i < d - 1; i++)
    for (j = 1; j <= d; j++)
      for (k = 0; k < deg; k++)
	if ((i+1)*d + j + k<= sink - 1)
	  AddEdge(i*d + j, (i+1)*d + j + k, RandomInteger(1, r), G);


  return G;
}


Graph *BasicLine(n, m, deg)
int n, m, deg;
{
  Graph *G;
  int i, j, source, sink, *x;

  x = (int *) myAlloc(deg * sizeof(int));

  G = (Graph *) myAlloc(sizeof(Graph));

  InitGraph(G, n*m + 2);  

  for (i = 0; i <= n*m + 1; i++)
    AddVertex(i, G);

  source = 0;
  sink = n*m + 1;

  for (i = 1; i <= m; i++){
    AddEdge(source, source + i, deg*MAX_RANGE, G);
    AddEdge(sink - i, sink, deg*MAX_RANGE, G);
  }

  for (i = source + 1; i < sink; i++){
      RandomSubset(1, m*deg, deg, x);
      for (j = 0; j < deg; j++)
	if (i + x[j] < sink)
	  AddEdge(i, i + x[j], RandomInteger(1, MAX_RANGE), G);
    }

  return G;
}


/* I'm not sure if these next two are needed - a bug has been fixed, so
   that the ranges are valid.  I think the reason for the MAX_DEGREE
   constant was so that the index into the Range array would be valid.
   A new constant LOG_RANGE has been introduced - values out of range
   are truncated.
   */


Graph *ExponentialLine(n, m, deg)
int n, m, deg;
{
  Graph *G;
  int i, j, source, sink, *x, r;

  x = (int *) myAlloc(deg*sizeof(int));

  G = (Graph *) myAlloc(sizeof(Graph));

  InitGraph(G, n*m + 2);  

  for (i = 0; i <= n*m + 1; i++)
    AddVertex(i, G);

  source = 0;
  sink = n*m + 1;

  for (i = 1; i <= m; i++){
    AddEdge(source, source + i, deg*MAX_RANGE, G);
    AddEdge(sink - i, sink, deg*MAX_RANGE, G);
  }

  for (i = source + 1; i < sink; i++){
      RandomSubset(1, m*deg, deg, x);
      for (j = 0; j < deg; j++){
	r = (x[j] - 1) / m;
        if (r >= LOG_RANGE)          /* Crude bug fix */
          r = 0;
	if (i + x[j] < sink)
	  AddEdge(i, i + x[j], RandomInteger(1, Range[r]), G);
      }
    }

  return G;
}

Graph *DExponentialLine(n, m, deg)
int n, m, deg;
{
  Graph *G;
  int i, j, source, sink, *x, r;

  x = (int *) myAlloc(deg*sizeof(int));

  G = (Graph *) myAlloc(sizeof(Graph));

  InitGraph(G, n*m + 2);  

  for (i = 0; i <= n*m + 1; i++)
    AddVertex(i, G);

  source = 0;
  sink = n*m + 1;

  for (i = 1; i <= m; i++){
    AddEdge(source, source + i, deg*MAX_RANGE, G);
    AddEdge(sink - i, sink, deg*MAX_RANGE, G);
  }

  for (i = source + 1; i < sink; i++){
      RandomSubset(-m*deg, m*deg, deg, x);
      for (j = 0; j < deg; j++){
	r = Abs((x[j] - 1) / m);
        if (r >= LOG_RANGE)          /* Crude bug fix */
          r = 0;
	if (i + x[j] < sink && i + x[j] > source && x[j] != 0)
	  AddEdge(i, i + x[j], RandomInteger(1, Range[r]), G);
      }
    }

  return G;
}


Graph *DinicBadCase(n)
int n;
{
  Graph *G;
  int i;

  G = (Graph *) myAlloc(sizeof(Graph));

  InitGraph(G, n);  

  for (i = 0; i < n; i++)
    AddVertex(i, G);

  for (i = 0; i < n-1; i++){
    AddEdge(i, i+1, n, G);
  }

  for (i = 0; i < n-2; i++){
    AddEdge(i, n - 1, 1, G);
  }


  return G;
}

Graph *GoldBadCase(n)
int n;
{
  Graph *G;
  int i;

  G = (Graph *) myAlloc(sizeof(Graph));

  InitGraph(G, 3*n+3);  

  for (i = 0; i < 3*n+3; i++)
    AddVertex(i, G);

  AddEdge(0, 1, n, G);

  for (i = 2; i < n+2; i++){
    AddEdge(1, i, n, G);
    AddEdge(i, i+n, 1, G);
    AddEdge(i+n, 2*n+2, n, G);
  }

  for (i = 2*n+2; i < 3*n+2; i++)
    AddEdge(i, i+1, n, G);

  return G;
}


Graph *Cheryian(n, m, c)
int n, m, c;
{
  Graph *G;
  int i;

  G = (Graph *) myAlloc(sizeof(Graph));

  InitGraph(G, 4*m*c + 6 + n);  

  AddVertex(0, G);
  AddVertex(1, G);
  AddVertex(2, G);
  AddVertex(3, G);

  Gadget(0, 1, n, m, c, G);
  Gadget(0, 2, n, m, c, G);
  Gadget(1, 3, n, m, c, G);
  Gadget(2, 3, n, m, c, G);

  Bridge(1, 2, n, G);

  Sink(3, G);

  return G;
}

Gadget(a, b, n, m, c, G)
int a, b, n, m, c;
Graph *G;
{
  int i, j, v, w;

  v = b;
  for (i = 0; i < m; i++){
    for (j = 0; j < c; j++){
      w = v;
      v = NewVertex(G);
      AddEdge(v, w, VERY_BIG, G);
    }
    AddEdge(a, v, n, G);
  }
}

Bridge(a, b, n, G)
int a, b, n;
Graph *G;
{
  int i, v, w,  v1, v2;

  v1 = NewVertex(G);
  v2 = NewVertex(G);
  AddEdge(a, v1, n, G);
  AddEdge(v2, b, n, G);
  for (i = 0; i < n; i++){
    v = NewVertex(G);
    w = NewVertex(G);
    AddEdge(v1, v, n, G);
    AddEdge(w, v2, n, G);
    AddEdge(v, w, 1, G);
  }
}

Sink(k, G)
int k;
Graph *G;
{
  AddEdge(k, NewVertex(G), VERY_BIG, G);
}

NewVertex(G)
Graph *G;
{
  AddVertex(G->size, G);
  return G->size - 1;
}


/* manip.c */


InitGraph(G, n)
Graph *G;
int n;
{
  int i;

  G->max_n = n;
  G->V = (int *) myAlloc(n * sizeof(int));
  G->A = (Edge **) myAlloc(n * sizeof(Edge **));

  for (i = 0; i < n; i++){
    G->A[i] = (Edge *) 0;
    G->V[i] = FALSE;
  }
  G->size = 0;  
  G->max_v = -1;
}



AddVertex(v, G)
int v;
Graph *G;
{
  if (G->V[v] == TRUE)
    Barf("Vertex already present");

  G->V[v] = TRUE;
  G->size++;
  if (v > G->max_v)
    G->max_v = v;
}


AddEdge(v1, v2, a, G)
int v1, v2, a;
Graph *G;
{
  Edge *e1, *e2, *EdgeLookup();

  if (v1 == v2)
    Barf("No Loops");

  if ((e1 = EdgeLookup(v1, v2, G)) != (Edge *) 0){
    e1->c += a;
    return;
  }

  e1 = (Edge *) myAlloc(sizeof(Edge));
  e2 = (Edge *) myAlloc(sizeof(Edge));

  e1->mate = e2;
  e2->mate = e1;

  e1->next = G->A[v1];
  G->A[v1] = e1;
  e1->t = v1;
  e1->h = v2;
  e1->c = a;

  e2->next = G->A[v2];
  G->A[v2] = e2;
  e2->t = v2;
  e2->h = v1;
  e2->c = 0;
}

Edge *EdgeLookup(v1, v2, G)
int v1, v2;
Graph *G;
{
  Edge *e;

  e = G->A[v1];
  while (e != (Edge *) 0){
    if (e->h == v2)
      return e;
    e = e->next;
  }
  return (Edge *) 0;
}


/* Count the number of edges with positive capacity */
int EdgeCount(G)
Graph *G;
{
  int i, count;
  Edge *e;

  count = 0;
  for (i = 0; i <= G->max_v; i++){
    if (G->V[i] == FALSE)
      continue;
    e = G->A[i];
    while (e != (Edge *) 0){
      if (e->c > 0)
	count++;
      e = e->next;
    }
  }
  return count;
}
  



Barf(s)
char *s;
{
  fprintf(stderr, "%s\n", s);
  exit(-1);
}

   
int Abs(x)
int x;
{
  return (x > 0) ? x : -x;
}


/* random.c -- functions dealing with randomization.

	RandomPermutation
        RandomInteger
        InitRandom
	RandomSubset
*/




/* RandomPermutation -- contruct a random permutation of the array perm.  
   It is assumed that the length of perm is n.  The algorithm used makes
   a pass through the array, randomly switching elements of the array.
*/
RandomPermutation (perm, n)
int perm[], n;
{
    int i, j, t;
			
    for (i = 0; i < n - 1; i++){
        j = RandomInteger(i, n-1);	/* Swap the element perm[i] with
*/
	t = perm[i];			/* a random element from the range
*/
	perm[i] = perm[j];		/* i..n-1.
*/ 
	perm[j] = t;
    }
}

RandPerm (perm, n)
int perm[], n;
{
    int i, j, t;
	
    for (i = 0; i < n; i++)
        perm[i] = i;

    for (i = 0; i < n - 1; i++){
        j = RandomInteger(i, n-1);	/* Swap the element perm[i] with
*/
	t = perm[i];			/* a random element from the range
*/
	perm[i] = perm[j];		/* i..n-1.
*/ 
	perm[j] = t;
    }
}



/* RandomInteger -- return a random integer from the range low .. high.
*/
int RandomInteger (low, high)
int low, high;
{
    return random() % (high - low + 1) + low;
}

/* InitRandom -- If the seed is non-zero, the random number generator is
   initialized with the seed, giving a fixed sequence of "random" numbers.
   If the seed is zero, then the time of day is used to intialize the random

   number generator. 
*/
InitRandom (seed)
int seed;
{
    struct timeval tp;
   
    if (seed == 0){
        gettimeofday(&tp, 0);
        srandom(tp.tv_sec + tp.tv_usec);
    }
    else
	srandom(seed);
}

/* RandomSubset - return n distinct values, randomly selected between
high and low, the algorithm is inefficient if n is large - this could
be improved */
RandomSubset(low, high, n, x)
int low, high, n, *x;
{
  int i, j, r, flag;

  if (high - low + 1 < n)
    Barf("Invalid range for Random Subset");

  i = 0;
  while (i < n){
    r = RandomInteger(low, high);
    flag = 0;
    for (j = 0; j < i; j++)
      if (x[j] == r)
	flag = 1;
    if (flag == 0)
      x[i++] = r;
  }
  
}







void GraphOutput(G, f, s, t)
Graph *G;
int s, t;
FILE *f;
{
  int i;
  void WriteVertex();

  fprintf(f, "p max %d %d\n", G->size, EdgeCount(G));
  fprintf(f, "n %d s\n", s + 1);
  fprintf(f, "n %d t\n", t + 1);
  for (i = 0; i <= G->max_v; i++){
    WriteVertex(i, G, f);
  }

}






/* for file output */
void WriteVertex(v, G, f)
int v;
Graph *G;
FILE *f;
{
  Edge *e;

  e = G->A[v];
  while (e != (Edge *) 0){
    if (e->c > 0){
        fprintf(f, "a %d %d %d\n", e->t + 1, e->h + 1, e->c);
    }
    e = e->next;
  }
}

/* test for allocation failure */
void *myAlloc(n)
int n;
{
  void *ptr;

  ptr = (void *) malloc(n);
  if (ptr == (void *) 0)
    Barf("Memory exhausted");

  return ptr;



}









