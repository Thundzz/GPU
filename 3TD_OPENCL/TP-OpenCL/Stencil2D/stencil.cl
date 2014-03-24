/* Original version */
/*
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
*/

__kernel void stencil(__global float *in,
		      __global float *out)
{
  int x = get_global_id(0);
  int y = get_global_id(1);
  unsigned int xloc = get_local_id(0);
  unsigned int yloc = get_local_id(1);


  in += OFFSET;
  out += OFFSET;
  __local float tile[TILE+2][TILE+2];

  //Starting the copy of the tile we want to handle
  if(x == 0 || y== 0 || y  == LINESIZE || x== LINESIZE)
    {
    }
  else{
    unsigned int index_in = y *LINESIZE + x;
    tile[xloc+1][yloc+1] = in[index_in];
    if (yloc == 1 )
      tile[xloc][yloc-1] = (y -1)*LINESIZE  +x;
    if (yloc == TILE )
      tile[xloc][yloc+1] = (y +1)*LINESIZE  +x;
    if (xloc == 1 )
      tile[xloc-1][yloc] = y*LINESIZE  +x-1;
    if (xloc == TILE )
      tile[xloc+1][yloc] = y*LINESIZE  +x +1; 
  }

  // Waiting for the copy to end
  barrier(CLK_LOCAL_MEM_FENCE);
  

  // Starting the real treatment

  out[y * LINESIZE + x] =
    C0 * tile[xloc][yloc] +
    C1 * ( tile[xloc - 1][yloc] +
	   tile[xloc + 1][yloc] +
	   tile[xloc][yloc - 1] +
	   tile[xloc][yloc + 1] );
}
