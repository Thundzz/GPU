#include "atom.h"
#include "default_defines.h"
#include "device.h"
#include "ocl.h"
#include "seq.h"
#include "sotl.h"

#ifdef HAVE_LIBGL
#include "vbo.h"
#endif

#ifdef HAVE_OMP
#include "openmp.h"
#endif

#include <stdio.h>

static void release_kernels(sotl_device_t *dev)
{
    int i;

    for (i = 0; i < KERNEL_TAB_SIZE; ++i) {
        clReleaseKernel(dev->kernel[i]);
    }
}

void device_finalize(sotl_device_t *dev)
{
    int i;

    if (dev->compute != SOTL_COMPUTE_OCL) {
        /* Do not free OpenCL ressources when an other mode is used. */
        return;
    }

    for (i = 0; i < 2; ++i) {
        clReleaseMemObject(dev->pos_buffer[i]);
        clReleaseMemObject(dev->speed_buffer[i]);
    }

    clReleaseMemObject(dev->box_buffer);
    clReleaseMemObject(dev->calc_offset_buffer);
    clReleaseMemObject(dev->min_buffer);
    clReleaseMemObject(dev->max_buffer);
    clReleaseMemObject(dev->fake_min_buffer);
    clReleaseMemObject(dev->fake_max_buffer);
    clReleaseMemObject(dev->domain_buffer);

    /* Release memory allocated by kernel objects. */
    release_kernels(dev);

    clReleaseCommandQueue(dev->queue);
    clReleaseProgram(dev->program);
    clReleaseContext(dev->context);
}

#define ALLOC_BUF(buf, size, flags, name)                               \
    do {                                                                \
        cl_int err;                                                     \
        buf = clCreateBuffer(dev->context, flags, size, NULL, &err);    \
        check(err, "Failed to create "name" buffer.");                  \
        dev->mem_allocated += size;                                     \
    } while (0)

#define ALLOC_RO_BUF(buf, size,  name) \
    ALLOC_BUF(buf, size, CL_MEM_READ_ONLY, name)

#define ALLOC_RW_BUF(buf, size,  name) \
    ALLOC_BUF(buf, size, CL_MEM_READ_WRITE, name)

static void create_gl_buffers(sotl_device_t *dev)
{
    if (!dev->display) {
        /* Do not create GL buffers when this device is not used for display. */
        return;
    }

#ifdef HAVE_LIBGL
    cl_int err;

    /* Shared VBO buffer with OpenGL. */
    vbo_buffer = clCreateFromGLBuffer(dev->context, CL_MEM_READ_WRITE,
                                      vbovid, &err);
    check(err, "Failed to map VBO buffer.\n");

#ifdef _SPHERE_MODE_
    /* Create model buffer. */
    size_t size = 2 * vertices_per_atom * 3 * sizeof(float);
    ALLOC_RW_BUF(model_buffer, size, "model_buffer");
#endif
#endif
}

void device_create_buffers(sotl_device_t *dev)
{
    size_t size;

    /* Compute size of position and speed buffers. */
    size = atom_set_size(&dev->atom_set);

    /* Allocate memory for borders (only in multi devices). */
    size += atom_set_border_size(&dev->atom_set);   /* left */
    size += atom_set_border_size(&dev->atom_set);   /* right */

    /* Create position and speed buffers. */
    for (int i = 0; i < 2; ++i) {
        ALLOC_RW_BUF(dev->pos_buffer[i], size, "pos_buffer(i)");
        ALLOC_RW_BUF(dev->speed_buffer[i], size, "speed_buffer(i)");
    }

    /* Init current position and speed buffers. */
    dev->cur_pb = 0;
    dev->cur_sb = 0;

    /* Compute size of box buffers. */
    unsigned total_boxes = ALRND(2 * SCAN_WG_SIZE, dev->domain.total_boxes + 1);
    size = total_boxes * sizeof(int);

    /* Create box buffers. */
    ALLOC_RW_BUF(dev->box_buffer, size, "box_buffer");
    ALLOC_RW_BUF(dev->calc_offset_buffer, size, "calc_offset_buffer");

    /* Create min and max buffers. */
    size = 3 * sizeof(calc_t);
    ALLOC_RO_BUF(dev->min_buffer, size, "min_buffer");
    ALLOC_RO_BUF(dev->max_buffer, size, "max_buffer");
    ALLOC_RO_BUF(dev->fake_min_buffer, size, "fake_min_buffer");
    ALLOC_RO_BUF(dev->fake_max_buffer, size, "fake_max_buffer");

    /* Create domain buffer. */
    size = 4 * sizeof(int);
    ALLOC_RO_BUF(dev->domain_buffer, size, "domain_buffer");

    /* Create GL buffers for display. */
    create_gl_buffers(dev);
}

#define WRITE_BUF(buffer, cb, offset, ptr, name)                            \
    do {                                                                    \
        cl_int err;                                                         \
        err = clEnqueueWriteBuffer(dev->queue, buffer, CL_TRUE, offset, cb, \
                                   ptr, 0, NULL, NULL);                     \
        check(err, "Failed to write to "name" buffer.");                    \
    } while (0)

static void write_gl_buffers(sotl_device_t *dev)
{
#ifdef HAVE_LIBGL
    size_t cb;
#endif

    if (!dev->display) {
        /* Do not write GL buffers when this device is not used for display. */
        return;
    }

#ifdef HAVE_LIBGL
    ocl_acquire(dev);

    /* Write VBO buffer. */
    cb = nb_vertices * 3 * sizeof(float);
    WRITE_BUF(vbo_buffer, cb, 0, vbo_vertex, "vbo_buffer");

#ifdef _SPHERE_MODE_
    /* Write model buffer. */
    cb = 2 * vertices_per_atom * 3 * sizeof(float);
    WRITE_BUF(model_buffer, cb, 0, vertex_model, "model_buffer");
#endif

    ocl_release(dev);
#endif
}

void device_write_buffers(sotl_device_t *dev)
{
    size_t cb, size, size_border, offset;

    /* Compute the size and the offset in bytes of data to write. */
    cb   = sizeof(calc_t) * dev->atom_set.natoms;
    size = atom_set_size(&dev->atom_set) / 3;

    /* Get size of one border (ie. left) (only in multi devices). */
    size_border = atom_set_border_size(&dev->atom_set) / 3;

    /* Write positions. */
    offset = size_border;
    WRITE_BUF(*cur_pos_buf(dev), cb, offset, dev->atom_set.pos.x, "pos_buffer(x)");
    offset += size + size_border * 2; /* for left and right borders. */
    WRITE_BUF(*cur_pos_buf(dev), cb, offset, dev->atom_set.pos.y, "pos_buffer(y)");
    offset += size + size_border * 2;
    WRITE_BUF(*cur_pos_buf(dev), cb, offset, dev->atom_set.pos.z, "pos_buffer(z)");

    /* Write back speeds. */
    offset = size_border;
    WRITE_BUF(*cur_spd_buf(dev), cb, offset, dev->atom_set.speed.dx, "speed_buffer(x)");
    offset += size + size_border * 2;
    WRITE_BUF(*cur_spd_buf(dev), cb, offset, dev->atom_set.speed.dy, "speed_buffer(y)");
    offset += size + size_border * 2;
    WRITE_BUF(*cur_spd_buf(dev), cb, offset, dev->atom_set.speed.dz, "speed_buffer(z)");

    /* Write min and max buffers. */
    cb = 3 * sizeof(calc_t);
    WRITE_BUF(dev->min_buffer, cb, 0, dev->domain.min_ext, "min_buffer");
    WRITE_BUF(dev->max_buffer, cb, 0, dev->domain.max_ext, "max_buffer");
    WRITE_BUF(dev->fake_min_buffer, cb, 0, dev->domain.min_border, "fake_min_buffer");
    WRITE_BUF(dev->fake_max_buffer, cb, 0, dev->domain.max_border, "fake_max_buffer");

    /* Write domain buffer. */
    cb = 4 * sizeof(int);
    WRITE_BUF(dev->domain_buffer, cb, 0, dev->domain.boxes, "domain_buffer");

    /* Write GL buffers for display. */
    write_gl_buffers(dev);
}

#define READ_BUF(buffer, cb, offset, ptr, name)                             \
    do {                                                                    \
        cl_int err;                                                         \
        err = clEnqueueReadBuffer(dev->queue, buffer, CL_TRUE, offset, cb,  \
                                  ptr, 0, NULL, NULL);                      \
        check(err, "Failed to read "name" back to host memory.");           \
        clFinish(dev->queue);                                               \
    } while (0)

void device_read_buffers(sotl_device_t *dev)
{
    device_read_back_pos(dev, dev->atom_set.pos.x, dev->atom_set.pos.y,
                         dev->atom_set.pos.z);
    device_read_back_spd(dev, dev->atom_set.speed.dx, dev->atom_set.speed.dy,
                         dev->atom_set.speed.dz);
}

void device_read_back_pos(sotl_device_t *dev, calc_t *pos_x, calc_t *pos_y,
                          calc_t *pos_z)
{
    size_t cb, size, size_border, offset;

    /* Compute the size and the offset in bytes of data to read. */
    cb   = sizeof(calc_t) * dev->atom_set.natoms;
    size = atom_set_size(&dev->atom_set) / 3;

    /* Get size of one border (ie. left) (only in multi devices). */
    size_border = atom_set_border_size(&dev->atom_set) / 3;

    /* Read back positions. */
    offset = size_border;
    READ_BUF(*cur_pos_buf(dev), cb, offset, pos_x, "pos_buffer(x)");
    offset += size + size_border * 2;
    READ_BUF(*cur_pos_buf(dev), cb, offset, pos_y, "pos_buffer(y)");
    offset += size + size_border * 2;
    READ_BUF(*cur_pos_buf(dev), cb, offset, pos_z, "pos_buffer(z)");
}

void device_read_back_spd(sotl_device_t *dev, calc_t *spd_x, calc_t *spd_y,
                          calc_t *spd_z)
{
    size_t cb, size, size_border, offset;

    /* Compute the size and the offset in bytes of data to read. */
    cb   = sizeof(calc_t) * dev->atom_set.natoms;
    size = atom_set_size(&dev->atom_set) / 3;

    /* Get size of one border (ie. left) (only in multi devices). */
    size_border = atom_set_border_size(&dev->atom_set) / 3;

    /* Read back speeds. */
    offset = size_border;
    READ_BUF(*cur_spd_buf(dev), cb, offset, spd_x, "speed_buffer(x)");
    offset += size + size_border * 2;
    READ_BUF(*cur_spd_buf(dev), cb, offset, spd_y, "speed_buffer(y)");
    offset += size + size_border * 2;
    READ_BUF(*cur_spd_buf(dev), cb, offset, spd_z, "speed_buffer(z)");
}

void device_one_step_move(sotl_device_t *dev)
{
    switch (dev->compute) {
        case SOTL_COMPUTE_OCL:
            ocl_one_step_move(dev);
            break;
        case SOTL_COMPUTE_SEQ:
            seq_one_step_move(dev);
            break;
#ifdef HAVE_OMP
        case SOTL_COMPUTE_OMP:
            omp_one_step_move(dev);
            break;
#endif
        default:
            sotl_log(ERROR, "Undefined compute method %d.\n", dev->compute);
            break;
    }
}

void cl_create_kernels(sotl_device_t *dev)
{
    cl_int err;
    int i;

    for (i = 0; i < KERNEL_TAB_SIZE; ++i) {
        dev->kernel[i] = clCreateKernel(dev->program, kernel_name(i), &err);
        check(err, "Failed to create a kernel object for '%s'.\n", kernel_name(i));
    }
}

void device_query_kernels_info(sotl_device_t *dev)
{
    cl_int err;
    int i;

    for (i = 0; i < KERNEL_TAB_SIZE; ++i) {
        sotl_log(INFO, "Work group info for kernel '%s'.\n", kernel_name(i));

        /* Work group size. */
        size_t work_group_size;
        err = clGetKernelWorkGroupInfo(dev->kernel[i], dev->id,
                                       CL_KERNEL_WORK_GROUP_SIZE,
                                       sizeof(size_t), &work_group_size, NULL);
        check(err, "Failed to get kernel work group size.\n");
        sotl_log(DEBUG, "WORK_GROUP_SIZE = %d\n", work_group_size);

        /* Compile work group size. */
        size_t compile_wgs[3];
        err = clGetKernelWorkGroupInfo(dev->kernel[i], dev->id,
                                       CL_KERNEL_COMPILE_WORK_GROUP_SIZE,
                                       sizeof(size_t) * 3, compile_wgs, NULL);
        check(err, "Failed to get kernel compile work group size.\n");
        sotl_log(DEBUG, "COMPILE_WORK_GROUP_SIZE = [%d, %d, %d]\n",
                 compile_wgs[0], compile_wgs[1], compile_wgs[2]);

        /* Local mem size. */
        cl_ulong local_mem_size;
        err = clGetKernelWorkGroupInfo(dev->kernel[i], dev->id,
                                       CL_KERNEL_LOCAL_MEM_SIZE,
                                       sizeof(cl_ulong), &local_mem_size, NULL);
        check(err, "Failed to get kernel work group size.\n");
        sotl_log(DEBUG, "LOCAL_MEM_SIZE = %d\n", local_mem_size);

        /* Preferred work group size multiple. */
        size_t preferred_wgs_multiple;
        err = clGetKernelWorkGroupInfo(dev->kernel[i], dev->id,
                                       CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
                                       sizeof(size_t), &preferred_wgs_multiple, NULL);
        check(err, "Failed to get kernel preferred work group size multiple.\n");
        sotl_log(DEBUG, "PREFERRED_WORK_GROUP_SIZE_MULTIPLE = %d\n", preferred_wgs_multiple);

        /* Private mem size. */
        cl_ulong private_mem_size;
        err = clGetKernelWorkGroupInfo(dev->kernel[i], dev->id,
                                       CL_KERNEL_PRIVATE_MEM_SIZE,
                                       sizeof(cl_ulong), &private_mem_size, NULL);
        check(err, "Failed to get kernel private mem size.\n");
        sotl_log(DEBUG, "PRIVATE_MEM_SIZE = %d\n", private_mem_size);
    }
}

static unsigned read_natoms_from_boxes(const sotl_device_t *dev, const int z)
{
    const int shift_z = dev->domain.boxes[0] * dev->domain.boxes[1];
    const size_t cb = sizeof(int);
    unsigned natoms = 0;
    size_t offset;

    if (z > 1) {
        offset = shift_z * z * cb;
        READ_BUF(dev->box_buffer, cb, offset, &natoms, "box_buffer[z - 1]");
    }

    return natoms;
}

static unsigned get_natoms_in_z(const sotl_device_t *dev, const int begin_z,
                                const int end_z)
{
    return read_natoms_from_boxes(dev, end_z)
           - read_natoms_from_boxes(dev, begin_z);
}

unsigned device_get_natoms_left(const sotl_device_t *dev)
{
    /* XXX: This function can be called just after scan to improve performance. */
    clFinish(dev->queue);

    /* Get number of atoms for the first two Z planes. */
    return get_natoms_in_z(dev, 0, 2);
}

unsigned device_get_natoms_right(const sotl_device_t *dev)
{
    /* XXX: This function can be called just after scan to improve performance. */
    clFinish(dev->queue);

    /* Get number of atoms for the last two Z planes. */
    return get_natoms_in_z(dev, dev->domain.boxes[2] - 2, dev->domain.boxes[2]);
}
