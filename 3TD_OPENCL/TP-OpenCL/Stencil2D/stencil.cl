
__kernel void stencil(__global float *in,
		      __global float *out)
{
   int x = get_global_id(0);
   int y = get_global_id(1);

   in += OFFSET;
   out += OFFSET;

   out[y * LINESIZE + x] =
     C0 * in[y * LINESIZE + x ] +
     C1 * ( in[y * LINESIZE + x - 1] +
	    in[y * LINESIZE + x + 1] +
	    in[(y - 1) * LINESIZE + x] +
	    in[(y + 1) * LINESIZE + x] );
}
