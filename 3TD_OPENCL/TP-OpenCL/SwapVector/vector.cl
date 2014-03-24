

__kernel void swap(__global float *vec,
			__global float *res)
{			
  __local float tmp[TILESIZE];

  int index = get_global_id(0);
  int nb = get_global_size(0);

  tmp[TILESIZE - get_local_id(0) -1] = vec[index];
  int grp = get_group_id(0);
  int nb_group = get_num_groups(0);
  
  barrier(CLK_LOCAL_MEM_FENCE);

  res[(nb_group - grp -1)*TILESIZE + get_global_id(0)] = tmp[get_local_id(0)];

}
