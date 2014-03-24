/* Original version*/

/*
__kernel void heat(__global float *next,
     __global float *previous)
{
   unsigned int i = get_global_id(0);
   
   next[i+1] = (previous[i] + previous[i+1] * 2 + previous[i+2]) / 4;
}
*/

__kernel void heat(__global float *next,
     __global float *previous)
{

  __local float tile [TILE +2];
   unsigned int i = get_global_id(0);
   
   unsigned int xloc = get_local_id(0);
   tile[xloc] = previous[i];
   if(xloc >= TILE -2)
     tile[xloc+2] = previous[i+2];

   
  // Waiting for every thread to complete his job
  barrier(CLK_LOCAL_MEM_FENCE);

   next[i+1] = (tile[xloc] + tile[xloc+1] * 2 + tile[xloc+2]) / 4;
}

