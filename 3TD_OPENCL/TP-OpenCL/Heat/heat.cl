
__kernel void
heat(__global float *next,
     __global float *previous)
{
   unsigned int i = get_global_id(0);

   next[i+1] = (previous[i] + previous[i+1] * 2 + previous[i+2]) / 4;
}
