#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <sys/time.h>

#include "tsp-types.h"
#include "tsp-job.h"


/* macro de mesure de temps, retourne une valeur en �secondes */
#define TIME_DIFF(t1, t2) \
        ((t2.tv_sec - t1.tv_sec) * 1000000 + (t2.tv_usec - t1.tv_usec))



/* variables globales */

/* arguments du programme */
int    Argc ;
char   **Argv ;

/* dernier minimum trouv� */
int    minimum ;

/* liste des t�ches */
TSPqueue     q ;

/* tableau des distances */
DTab_t    distance ;

/* nombre de villes */
int NrTowns ;

#define MAXX	100
#define MAXY	100
typedef struct
		{
		 int x, y ;
		} coor_t ;

typedef coor_t coortab_t [MAXE] ;
coortab_t towns ;

#define NUM 3

/* initialisation du tableau des distances */
void genmap () 
{
int i, j ;
int dx, dy;

 NrTowns = atoi (Argv[1]) ;
 if (NrTowns > MAXE) {
  fprintf(stderr,"trop de villes, augmentez MAXE dans tsp-types.h");
  exit(1);
 }

 srand (atoi(Argv[2])) ;

 for (i=0; i<NrTowns; i++)
  {
   towns [i].x = rand () % MAXX ;
   towns [i].y = rand () % MAXY ;
  }

 for (i=0; i<NrTowns; i++)
  {
   for (j=0; j<NrTowns; j++)
    {
     /* Un peu r�aliste */
     dx = towns [i].x - towns [j].x ;
     dy = towns [i].y - towns [j].y ;
     distance [i][j] = (int) sqrt ((double) ((dx * dx) + (dy * dy))) ;
    }
  }
}
 

/* impression tableau des distances, pour v�rifier au besoin */
void PrintDistTab ()
{
 int i, j ;

 printf ("NrTowns = %d\n", NrTowns) ;

 for (i=0; i<NrTowns; i++)
  {
   printf ("distance [%1d]",i) ;
   for (j=0; j<NrTowns; j++)
    {
     printf (" [%2d:%2d] ", j, distance[i][j]) ;
    }
   printf (";\n\n") ;
  }
 printf ("done ...\n") ;

}

void printPath (Path_t path)
{
 char toprint[MAXY][MAXX+1];
 int i;
 int x, y;

 memset(toprint, ' ', sizeof(toprint));
 for (i = 0; i < NrTowns-1; i++)
 {
  int x1 = towns[path[i]].x;
  int y1 = towns[path[i]].y;
  int x2 = towns[path[i+1]].x;
  int y2 = towns[path[i+1]].y;

  if (abs(x2-x1) > abs(y2-y1))
   for (x = 1; x < abs(x2-x1); x++) {
    toprint[y1+x*(y2-y1)/abs(x2-x1)][x1+x*(x2-x1)/abs(x2-x1)] = '*';
   }
  else
   for (y = 1; y < abs(y2-y1); y++) {
    toprint[y1+y*(y2-y1)/abs(y2-y1)][x1+y*(x2-x1)/abs(y2-y1)] = '*';
   }
 }

 for (i = 0; i < NrTowns; i++)
  toprint[towns[i].y][towns[i].x] = '#';
 for (y = 0; y < MAXY; y++) {
  toprint[y][MAXX] = 0;
  printf("%s\n", toprint[y]);
 }
}

/* r�solution du probl�me du voyageur de commerce */

int present (int city, int hops, Path_t path)
{
 unsigned int i ;

 for (i=0; i<hops; i++)
   if (path [i] == city) return 1 ;
 return 0 ;
}

void tsp (int hops, int len, Path_t path, int *cuts)
{
 int i ;
 int me, dist ;

 if (len >= minimum)
   {
     cuts[hops]++ ;
     return;
   }

 if (hops == NrTowns)
   {
     if (len + distance[0][path[NrTowns-1]] < minimum)
      {
       minimum = len + distance[0][path[NrTowns-1]] ;
       fprintf (stderr, "found path len = %3d :", minimum) ;
       for (i=0; i < NrTowns; i++)
         fprintf (stderr, "%2d ", path[i]) ;
       fprintf (stderr, "\n") ;
       //printPath(path);
      }
   }
 else
   {
   me = path [hops-1] ;

   for (i=0; i < NrTowns; i++)
     {
     if (!present (i, hops, path))
       {
         path [hops] = i ;
         dist = distance[me][i] ;
         tsp (hops+1, len+dist, path, cuts) ;
       }
     }
 
  }
}

void tsp_partiel (int hops, int len, Path_t path, int *cuts)
{
 int i ;
 int me, dist ;

 if (len >= minimum)
   {
     cuts[hops]++ ;
     return;
   }

 if (hops == NUM)
   {
     add_job (&q, path, hops, len);
   }
 else
   {
   me = path [hops-1] ;

   for (i=0; i < NrTowns; i++)
     {
     if (!present (i, hops, path))
       {
         path [hops] = i ;
         dist = distance[me][i] ;
         tsp_partiel (hops+1, len+dist, path, cuts) ;
       }
     }
 
  }
}

int main (int argc, char **argv)
{
   unsigned long temps;
   int i;
   struct timeval t1, t2;

   int *cuts = NULL; /* Juste � des fins de statistiques pour savoir combien de fois on a pu optimiser */

   if (argc != 3)
     {
	fprintf (stderr, "Usage: %s  <ncities> <seed> \n",argv[0]) ;
	exit (1) ;
     }
   
   Argc = argc ;
   Argv = argv ;
   minimum = INT_MAX ;
      
   printf ("ncities = %3d\n", atoi(argv[1])) ;
   
   init_queue (&q);
   genmap () ; 

   gettimeofday(&t1,NULL);

   {
      Path_t path;
      int i;

      for(i=0;i<MAXE;i++)
        path[i] = -1 ;
      path [0] = 0;
      
      tsp_partiel (1, 0, path, cuts);

      no_more_jobs (&q);
   }
    
   gettimeofday(&t2,NULL);

   temps = TIME_DIFF(t1,t2);
  
   printf("time = %ld.%03ldms\n", temps/1000, temps%1000);
   {
      Path_t path;
      int hops, len, n = 0;

      cuts = calloc(NrTowns+1,sizeof(*cuts));

      while (!empty_queue(&q))
	 {
	    n++;
	    get_job(&q, path, &hops, &len);
	    tsp(hops, len, path, cuts);
	 }
      printf("%d jobs", n);
   }
    
   gettimeofday(&t2,NULL);

   temps = TIME_DIFF(t1,t2);
   printf("time = %ld.%03ldms (cuts :", temps/1000, temps%1000);
   for (i=0; i<NrTowns; i++)
     printf(" %d",cuts[i]);
   printf(")\n");

   return 0 ;
}
