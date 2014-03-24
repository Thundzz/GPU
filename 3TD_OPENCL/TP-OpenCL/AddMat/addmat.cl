
__kernel void addmat(__global float *A,
		     __global float *B,
		     __global float *C)
{
  const int x = get_global_id(0);
  const int y = get_global_id(1);
  const int offset = get_global_size(0) * y + x;

  C[offset] = A[offset] + B[offset];
}
