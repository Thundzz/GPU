
__kernel void vecmul(__global float *vec,
		     __global float *res,
		     float k)
{
  int index = get_global_id(0);

   res[index] = vec[index] * k;
}
