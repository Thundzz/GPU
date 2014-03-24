/* Original naive kernel */

/*
__kernel void transpose(__global float *in,
			__global float *out)
{
  int x = get_global_id(0);
  int y = get_global_id(1);
  __local float * tile;
		     
  out[x * get_global_size(0) + y] = in[y * get_global_size(0) + x];
}

*/


__kernel void transpose(__global float *in,
			__global float *out)
{
  int x = get_global_id(0);
  int y = get_global_id(1);
  __local float tile[(TILE+1)*(TILE+1)];
 

 //Transposing into the local tile buffer  
  
  unsigned int index_in = y *SIZE  + x;
  unsigned int local_index = get_local_id(1)*(TILE)+get_local_id(0);

  tile[local_index] = in[index_in];



		     
  // Waiting for every thread to complete his job 
  barrier(CLK_LOCAL_MEM_FENCE);

  // Copying the content of tile local tile into the output memory buffer

  
  unsigned int xIndex = get_group_id(1) * TILE + get_local_id(0);
  unsigned int yIndex = get_group_id(0) * TILE + get_local_id(1);

  unsigned int inverted_local_index = get_local_id(0)*(TILE)+get_local_id(1);
  unsigned int index_out = yIndex * SIZE + xIndex;
  out[index_out] = tile[inverted_local_index];

}
