#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <stdint.h>

int N=2048;
int P=8;
int dyn=0;


pthread_mutex_t dynaMutex = PTHREAD_MUTEX_INITIALIZER;

typedef union u_tick {
        unsigned long long tick;
        struct {
                unsigned low;
                unsigned high;
        };
} tick_t;


#define GET_TICK(t) \
           __asm__ volatile("rdtsc" : "=a" ((t).low), "=d" ((t).high))
#define TICK_DIFF(t1, t2) \
           ((t2).tick - (t1).tick)



int (*fun)(int);

int fcroissant (int i)
{

  volatile int s = 0,x;
  int max = i/4;
  for(x=0; x < max ; x++)
    {
      int y;
      for(y=0; y < max ; y++)
	s+=y;
    }
  return s;
}

int fpetit (int i)
{
  return 1;
}


int fconstant (int i)
{
  volatile int s = 0,x;
  int max = 50;
  for(x=0; x < max ; x++)
    {
      int y;
      for(y=0; y < max ; y++)
	s+=y;
    }
  return s;
}

int fperiodique (int i)
{
  volatile int s = 0,x;
  int max = i%200;
  if (max<2)
    max=5000;

  for(x=0; x < max ; x++)
    {
      int y;
      for(y=0; y < max ; y++)
	s+=y;
    }
  return s;
}


long long sequentiel()
{
  int i;
  tick_t t1,t2;

  GET_TICK(t1);

  for(i=0; i< N; i++)
    fun(i);

  GET_TICK(t2);
  
  return TICK_DIFF(t1,t2);
}

void *cyclique(void *param)
{
	int * hax = (int *) param;
	int threadNum = hax[0];

  for(int i=threadNum; i<N; i+=P)
    fun(i);

  return NULL;
}

void *bloc(void *param)
{
	int * hax = (int *) param;
	int threadNum = hax[0];

  for(int i=threadNum*(N/P); i<N && i<(threadNum+1)*(N/P); i++)
    fun(i);

  return NULL;
}
// __sync_Fetch_and_add(&index, 1);
void *dynamique(void *param)
{
	int maTache;
	while (dyn< N)
	{
    pthread_mutex_lock(&dynaMutex);
		maTache = dyn ++ ;
		pthread_mutex_unlock(&dynaMutex);
		fun (maTache);
			
		//fprintf(stderr, "%d\n", maTache);
	}

  return NULL;
}


long long lance(void *(thread_fun)(void *param))
{
	dyn = 0;
  tick_t t1,t2;
  pthread_t t[P];
  int id[P];
  //int returnValues[P];
  for(int i = 0 ; i < P ; i++)
  {
		id[i] = i;
	}

  
  GET_TICK(t1);
  
  
  for (int i = 0; i<P; i++){
		pthread_create(&t[i], NULL, thread_fun, &id[i]);
		//fprintf(stderr, "%d\n",id[i]);
	}
  
  for (int i = 0; i<P; i++){
		pthread_join(t[i], NULL /* &returnValues[i]*/);
  }
  
  GET_TICK(t2);
  
  return TICK_DIFF(t1,t2);
}

void benchmark(char *s)
{
  float f, b,c,d;

  printf("%s (P=%d, N=%d)\n",s, P, N);

  f = ((float)sequentiel())/1000000;
  printf("sequentiel: \t%.3f Mcycles\n",f);

  b = ((float)lance(bloc))/1000000;
  printf("bloc:   \t%.3f Mcycles\t accélération %.1f\n",b,f/b);
  
  c = ((float)lance(cyclique))/1000000;
  printf("cyclique: \t%.3f Mcycles\t accélération %.1f\n",c,f/c);
  
  d = ((float)lance(dynamique))/1000000;
  printf("dynamique: \t%.3f Mcycles\t accélération %.1f\n",d,f/d);
}

int main(int argc, char **argv){	

  if(argc != 2)
    {
      fprintf(stderr,"%s nbthreads\n",argv[0]);
      exit(1);
    }
  
  P=atoi(argv[1]);

  fun = fpetit;
  benchmark("    PETIT");

  fun = fconstant;
  benchmark("    CONSTANT");

  fun = fcroissant;
  benchmark("    CROISSANT"); 
  
  fun = fperiodique;
  benchmark("    PERIODIQUE"); 
  
  return 0;
}
