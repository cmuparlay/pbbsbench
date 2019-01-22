/* generator of hard maxflow problems */
/* 01/09/94 - Stanford Computer Science Department */
/* Boris Cherkassky - cher@theory.stanford.edu, on.cher@zib-berlin.de */
/* Andrew V. Goldberg - goldberg@cs.stanford.edu */

#include <stdio.h>
#include <stdlib.h>

main ( argc, argv )

int argc;
char* argv[];

{

int n, i, d;

if ( argc < 2 ) goto usage;

n = atoi(argv[1]);

if ( n < 2 ) goto usage;

printf("c very bad maxflow problem\n");

printf("p max %d %d\n", 4*n+6, 6*n+7);

printf("n 1 s\n");
printf("n 2 t\n");

/* first terrible graph */

for ( i = 0; i<n ; i++ )
   {
     printf("a %d %d %d\n", i+3, i+4, n-i+1);
     printf("a %d %d %d\n", i+3, n+4, 1);
   }

printf("a %d %d %d\n", n+3, 2*n+4, 1);
printf("a %d %d %d\n", n+3, n+4, 1);

for ( i = n+3; i <= 2*n+2; i ++ )
   printf("a %d %d %d\n", i+1, i+2, n+1);

/* second horrible graph */

d = 2*n+4;

for ( i = d; i<=2*n+d ; i++ )
     printf("a %d %d %d\n", i+1, i+2, n);

for ( i = 0; i < n; i ++ )
     printf("a %d %d %d\n", i+d+1, 2*n+2-i+d, 1);


/* edges from source and to sink */

printf("a 1 3 1000000\n");
printf("a 1 %d 1000000\n", d+1);
printf("a %d 2 1000000\n", d);
printf("a %d 2 1000000\n", 4*n+6);


exit(0);

 usage:

printf(" usage: ak n (n>1)\n");
}
