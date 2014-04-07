
#include "default_defines.h"
#include "kernel_list.h"

static char *kernels_name[KERNEL_TAB_SIZE] = {
    "gravity",
    "eating",
    "growing",
    "reset_int_buffer",
    "box_count_all_atoms",
    "box_count_own_atoms",
    "scan",
    "scan_down_step",
    "copy_buffer",
    "box_sort_all_atoms",
    "box_sort_own_atoms",
    "I_believe_I_can_slide",
#ifdef MACHINE_DE_PAUVRE
    "basic_lennard_jones",
#else
    "n2_lennard_jones",
#endif
    "border_collision", 
    "update_position",
#ifdef _SPHERE_MODE_
    "update_vertices",
#else
    "update_vertice",
#endif

    "zero_speed",
    "atom_collision",
    "null_kernel",
    "reset_calc_t_buffer",
};



char *kernel_name(unsigned kernel_num)
{
  return kernels_name[kernel_num];
}
