#include "default_defines.h"
#include "global_definitions.h"
#include "device.h"
#include "seq.h"
#include "sotl.h"

#ifdef HAVE_LIBGL
#include "vbo.h"
#endif

#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <time.h>

static int *atom_state = NULL;

static int n = 0;
static float time_force = 0.0;
static float time_bounce = 0.0;
static float time_move = 0.0;
static float time_gravity = 0.0;

#ifdef HAVE_LIBGL

#define SHOCK_PERIOD  50
#define TIME_DIFF(t1, t2) \
		((t2.tv_sec - t1.tv_sec) * 1000000 + (t2.tv_usec - t1.tv_usec))

// Update OpenGL Vertex Buffer Object
//
static void seq_update_vbo (sotl_device_t *dev)
{
    sotl_atom_set_t *set = &dev->atom_set;
    //sotl_domain_t *domain = &dev->domain;

    for (unsigned n = 0; n < set->natoms; n++)
    {
        vbo_vertex[n*3 + 0] = set->pos.x[n];
        vbo_vertex[n*3 + 1] = set->pos.y[n];
        vbo_vertex[n*3 + 2] = set->pos.z[n];


        // Atom color depends on z coordinate
        if(atom_state[n]!= 0)
        {
            float ratio =(float) ( atom_state[n]) / (float) SHOCK_PERIOD;

            vbo_color[n*3 + 0] = (1.0 - ratio) * atom_color[0].R + ratio * 1.0;
            vbo_color[n*3 + 1] = (1.0 - ratio) * atom_color[0].G + ratio * 0.0;
            vbo_color[n*3 + 2] = (1.0 - ratio) * atom_color[0].B + ratio * 0.0;
            atom_state[n]--;
        }
    }
}
#endif

// Update positions of atoms by adding (dx, dy, dz)
//
static void seq_move (sotl_device_t *dev)
{
    struct timeval tv1,tv2;
    sotl_atom_set_t *set = &dev->atom_set;
    gettimeofday(&tv1,NULL);

    for (unsigned n = 0; n < set->natoms; n++)
    {
        set->pos.x[n] += set->speed.dx[n];
        set->pos.y[n] += set->speed.dy[n];
        set->pos.z[n] += set->speed.dz[n];
    }
    gettimeofday(&tv2,NULL);
    time_move += (float)TIME_DIFF(tv1,tv2);
}

// Apply gravity force
//
static void seq_gravity (sotl_device_t *dev)
{
    struct timeval tv1,tv2;
    gettimeofday(&tv1,NULL);
    sotl_atom_set_t *set = &dev->atom_set;
    const calc_t g = 0.005;

    //TODO
    for (unsigned n = 0; n < set->natoms; n++)
    {
        set->speed.dz[n] -= g;
    }
    gettimeofday(&tv2,NULL);
    time_gravity += (float)TIME_DIFF(tv1,tv2);
}

static void seq_bounce (sotl_device_t *dev)
{
    struct timeval tv1,tv2;
    // EN COURS
    sotl_atom_set_t *set = &dev->atom_set;
    //  sotl_domain_t *domain = &dev->domain;
    float xmin = dev->domain.min_ext[0];
    float ymin = dev->domain.min_ext[1];
    float zmin = dev->domain.min_ext[2];

    float xmax = dev->domain.max_ext[0];
    float ymax = dev->domain.max_ext[1];
    float zmax = dev->domain.max_ext[2];

    gettimeofday(&tv1,NULL);
    for (unsigned n = 0; n < set->natoms; n++)
    {

        if(set->pos.x[n] < xmin || set->pos.x[n] > xmax)
        {
            set->speed.dx[n] *= -1.0;
            atom_state[n]= SHOCK_PERIOD;
        }
        if(set->pos.y[n] < ymin || set->pos.y[n] > ymax)
        {
            set->speed.dy[n] *= -1.0;
            atom_state[n]= SHOCK_PERIOD;
        }
        if(set->pos.z[n] < zmin || set->pos.z[n] > zmax)
        {
            set->speed.dz[n] *= -1.0;
            atom_state[n]= SHOCK_PERIOD;
        }
    }
    gettimeofday(&tv2,NULL);
    time_bounce += (float)TIME_DIFF(tv1,tv2);
}


static calc_t squared_distance (sotl_atom_set_t *set, unsigned p1, unsigned p2)
{
    calc_t *pos1 = set->pos.x + p1,
            *pos2 = set->pos.x + p2;

    calc_t dx = pos2[0] - pos1[0],
           dy = pos2[set->offset] - pos1[set->offset],
           dz = pos2[set->offset*2] - pos1[set->offset*2];

    return dx * dx + dy * dy + dz * dz;
}

static calc_t lennard_jones (calc_t r2)
{
    calc_t rr2 = 1.0 / r2;
    calc_t r6;

    r6 = LENNARD_SIGMA * LENNARD_SIGMA * rr2;
    r6 = r6 * r6 * r6;

    return 24 * LENNARD_EPSILON * rr2 * (2.0f * r6 * r6 - r6);
}

static void seq_force (sotl_device_t *dev)
{
    struct timeval tv1,tv2;
    gettimeofday(&tv1,NULL);
    sotl_atom_set_t *set = &dev->atom_set;

    for (unsigned current = 0; current < set->natoms; current++)
    {
        calc_t force[3] = { 0.0, 0.0, 0.0 };

        for (unsigned other = 0; other < set->natoms; other++)
            if (current != other)
            {
                calc_t sq_dist = squared_distance (set, current, other);

                if (sq_dist < LENNARD_SQUARED_CUTOFF)
                {
                    calc_t intensity = lennard_jones (sq_dist);

                    force[0] += intensity * (set->pos.x[current] - set->pos.x[other]);
                    force[1] += intensity * (set->pos.x[set->offset + current] -
                                             set->pos.x[set->offset + other]);
                    force[2] += intensity * (set->pos.x[set->offset * 2 + current] -
                                             set->pos.x[set->offset * 2 + other]);
                }

            }

        set->speed.dx[current] += force[0];
        set->speed.dx[set->offset + current] += force[1];
        set->speed.dx[set->offset * 2 + current] += force[2];

    }
    gettimeofday(&tv2,NULL);
    time_force += (float)TIME_DIFF(tv1,tv2);
}


// Main simulation function
//
void seq_one_step_move (sotl_device_t *dev)
{
    // Apply gravity force
    //
    if (gravity_enabled)
        seq_gravity (dev);

    // Compute interactions between atoms
    //
    if (force_enabled)
        seq_force (dev);

    // Bounce on borders
    //
    if(borders_enabled)
        seq_bounce (dev);

    // Update positions
    //
    seq_move (dev);
    n++;
#ifdef HAVE_LIBGL
    // Update OpenGL position
    //
    if (dev->display)
        seq_update_vbo (dev);
#endif
}

void seq_init (sotl_device_t *dev)
{
#ifdef _SPHERE_MODE_
    sotl_log(ERROR, "Sequential implementation does currently not support SPHERE_MODE\n");
    exit (1);
#endif

    borders_enabled = 1;

    dev->compute = SOTL_COMPUTE_SEQ; // dummy op to avoid warning
}

void seq_alloc_buffers (sotl_device_t *dev)
{
    atom_state = calloc(dev->atom_set.natoms, sizeof(int));
    printf("natoms: %d\n", dev->atom_set.natoms);
}

void seq_finalize (sotl_device_t *dev)
{
    printf("Calculation time after %d iterations : %.1fms\n", n,(time_move+time_force+time_bounce)/1000);
    printf("Details per iteration : time_move %.1fus - time_force %.1fus - time_bounce %.1fus\n", time_move/n, time_force/n, time_bounce/n);
    free(atom_state);

    dev->compute = SOTL_COMPUTE_SEQ; // dummy op to avoid warning
}
