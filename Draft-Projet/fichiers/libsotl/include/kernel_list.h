#ifndef KERNEL_LIST_H
#define KERNEL_LIST_H


enum {
    KERNEL_GRAVITY,
    KERNEL_EATING_PACMAN,
    KERNEL_GROWING_GHOST,
    KERNEL_RESET_BOXES,
    KERNEL_BOX_COUNT_ALL_ATOMS,
    KERNEL_BOX_COUNT_OWN_ATOMS,
    KERNEL_SCAN,
    KERNEL_SCAN_DOWN_STEP,
    KERNEL_COPY_BUFFER,
    KERNEL_BOX_SORT_ALL_ATOMS,
    KERNEL_BOX_SORT_OWN_ATOMS,
    KERNEL_FORCE,
    KERNEL_FORCE_N2,
    KERNEL_BORDER, 
    KERNEL_UPDATE_POSTION,
    KERNEL_UPDATE_VERTICES,

    KERNEL_ZERO_SPEED,
    KERNEL_COLLISION, 
    KERNEL_NULL,
    KERNEL_RESET_CALC_T_BUFFER,

    KERNEL_TAB_SIZE
};

char *kernel_name(unsigned kernel_num);

#endif
