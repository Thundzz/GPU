#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define N 1024

#define TIME_DIFF(t1, t2) \
  ((double)((t2.tv_sec - t1.tv_sec) * 1000000 + (t2.tv_usec - t1.tv_usec)))

char a[N][N];
char b[N][N];
char c[N][N];
char d[N][N];

void mulMat1() {
	int i,j,k;
	for(j=0;j<N;j++) {
		for(i=0;i<N;i++) {
			c[i][j]=0;
			for(k=0;k<N;k++)
				c[i][j]+=a[i][k]*b[k][j];			
		}
	}
}

int compare()
{
	int true = 1;
	int i, j;
	for(i=0; i<N; i++)
		for(j=0; j< N; j++ )
			true = (true && (c[i][j]==d[i][j]));
	return true;
}

void mulMat2() {
	
	int i,j,k;
	for(i=0; i<N; i++)
		for(j=0; j< N; j++ )
			d[i][j]=0;
	for(i=0;i<N;i++) {	
		for(k=0;k<N;k++)
			for(j=0;j<N;j++) {
			{
				d[i][j]+=a[i][k]*b[k][j];			
			}
		}
	}	
}

int main(int argc, char ** argv){
        struct timeval t1,t2,t3;
	int i,j;
	for(i=0;i<N;i++)
		for(j=0;j<N;j++){			
			a[i][j]=rand()%2;
			b[i][j]=rand()%2;
		}

        gettimeofday(&t1,NULL);
	mulMat1();
        gettimeofday(&t2,NULL);
	mulMat2();
        gettimeofday(&t3,NULL);

	printf("%g / %g  = acceleration = %g\n", TIME_DIFF(t1,t2) / 1000,  TIME_DIFF(t2,t3) / 1000, TIME_DIFF(t1,t2)/TIME_DIFF(t2,t3));
	printf("is equal : %d\n", compare());
	return 0;
}
