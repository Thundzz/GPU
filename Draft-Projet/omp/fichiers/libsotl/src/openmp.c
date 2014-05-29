#include "default_defines.h"
#include "global_definitions.h"
#include "device.h"
#include "openmp.h"
#include "sotl.h"

#ifdef HAVE_LIBGL
#include "vbo.h"
#endif

#include <sys/time.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

int n = 0;
float time_force = 0.0;
float time_bounce = 0.0;
float time_move = 0.0;

static int *atom_state = NULL;

// Déclaration de variables pour le tri par boites.
static sotl_atom_set_t * tmp_atom_set;
static int * box_count;
static int * box_count_cummul;
static int * atom_box_index;

typedef struct index_list_t
{
    int elements[27];
    int nbElements;
    int current;
} index_list;

#ifdef HAVE_LIBGL

#define SHOCK_PERIOD  50
#define TIME_DIFF(t1, t2) \
		((t2.tv_sec - t1.tv_sec) * 1000000 + (t2.tv_usec - t1.tv_usec))

struct thread_color
{
    float R, G, B;
} proc_color[] =
{
    {1.0,0.0,0.0},
    {0.0,1.0,0.0},
    {0.0,0.0,1.0},
    {1.0,1.0,1.0},
    {1.0,1.0,0.0}
};


// Update OpenGL Vertex Buffer Object
//
static void omp_update_vbo (sotl_device_t *dev)
{
    sotl_atom_set_t *set = &dev->atom_set;
    //sotl_domain_t *domain = &dev->domain;
    #pragma omp parallel for
    for (unsigned n = 0; n < set->natoms; n++)
    {
        vbo_vertex[n*3 + 0] = set->pos.x[n];
        vbo_vertex[n*3 + 1] = set->pos.y[n];
        vbo_vertex[n*3 + 2] = set->pos.z[n];

        // Atom color depends on z coordinate
        /* if(atom_state[n]!= 0) */
        /* { */
        /*   float ratio =(float) ( atom_state[n]) / (float) SHOCK_PERIOD; */
        /*   vbo_color[n*3 + 0] = (1.0 - ratio) * atom_color[0].R + ratio * 1.0; */
        /*   vbo_color[n*3 + 1] = (1.0 - ratio) * atom_color[0].G + ratio * 0.0; */
        /*   vbo_color[n*3 + 2] = (1.0 - ratio) * atom_color[0].B + ratio * 0.0; */
        /*   atom_state[n]--; */
        /* } */
        /* Couleur en fonction du coeur
        {
            vbo_color[n*3 + 0] = proc_color[omp_get_thread_num() % 4].R;
            vbo_color[n*3 + 1] = proc_color[omp_get_thread_num() % 4].G;
            vbo_color[n*3 + 2] = proc_color[omp_get_thread_num() % 4].B;
        }*/
        /*Couleur en fonction de la boite
        {
            vbo_color[n*3 + 0] = proc_color[atom_box_index[n] % 5].R;
            vbo_color[n*3 + 1] = proc_color[atom_box_index[n] % 5].G;
            vbo_color[n*3 + 2] = proc_color[atom_box_index[n] % 5].B;
        }*/
    }
}
#endif


// Fonctions du tri par boites :


static void set_at_right_position(sotl_atom_set_t *orig, sotl_atom_set_t *dest, int orig_pos, int dest_pos)
{
    dest->pos.z[dest_pos] = orig->pos.z[orig_pos];
    dest->pos.y[dest_pos] = orig->pos.y[orig_pos];
    dest->pos.x[dest_pos] = orig->pos.x[orig_pos];
    dest->speed.dz[dest_pos] = orig->speed.dz[orig_pos];
    dest->speed.dy[dest_pos] = orig->speed.dy[orig_pos];
    dest->speed.dx[dest_pos] = orig->speed.dx[orig_pos];
}

void box_variables_init()
{
    if (tmp_atom_set != NULL)
        return;

    sotl_atom_set_t * glob_aset = get_global_atom_set();
    sotl_domain_t * glob_dom = get_global_domain();

    tmp_atom_set = (sotl_atom_set_t *) malloc(sizeof(sotl_atom_set_t));

    atom_set_init(tmp_atom_set, glob_aset->natoms, glob_aset->natoms);
    box_count = (int *) calloc (glob_dom->total_boxes , sizeof(int));
    box_count_cummul = (int *) calloc (glob_dom->total_boxes , sizeof(int));
    atom_box_index = (int *) calloc (glob_aset->natoms , sizeof(int));
}

void box_count_reinit()
{
    sotl_domain_t * glob_dom = get_global_domain();
    sotl_atom_set_t * glob_aset = get_global_atom_set();
    memset((void *) box_count, 0, (glob_dom->total_boxes)*sizeof(int));
    memset((void *) box_count_cummul, 0, (glob_dom->total_boxes)*sizeof(int));
    memset((void *) atom_box_index, 0, (glob_aset->natoms)*sizeof(int));
}

int atom_box_calc(int n)
{
    sotl_atom_set_t * glob_set = get_global_atom_set();
    sotl_domain_t * glob_dom = get_global_domain();

    int bx, by, bz;
    int box_id;

    bx = (int)((glob_set->pos.x[n] - glob_dom->min_border[0])*BOX_SIZE_INV);
    by = (int)((glob_set->pos.y[n] - glob_dom->min_border[1])*BOX_SIZE_INV);
    bz = (int)((glob_set->pos.z[n] - glob_dom->min_border[2])*BOX_SIZE_INV);
    /*
        if(bx <= 0)
            bx =0;
        if (bx >= glob_dom->boxes[0])
            bx = glob_dom->boxes[0]-1;

        if(by <= 0)
            by =0;
        if (by >= glob_dom->boxes[1])
            by = glob_dom->boxes[1]-1;

        if(bz <= 0)
            bz =0;
        if (bz >= glob_dom->boxes[2])
            bz = glob_dom->boxes[2]-1;
    */
    box_id =  bz *(glob_dom->boxes[0] * glob_dom->boxes[1]) +
              by * glob_dom->boxes[0] +
              bx;
    fprintf(stderr, "%d %g %g\n", bx,glob_set->pos.x[n], glob_dom->min_border[0]);
    assert( (unsigned)box_id <= glob_dom->total_boxes);
    return (int)box_id;
}


void atom_set_box_calc()
{
    int box_id;
    sotl_atom_set_t * glob_set = get_global_atom_set();

    for (unsigned n = 0; n < glob_set->natoms; n++)
    {
        box_id =  atom_box_calc(n);
        atom_box_index[n] = box_id;
        box_count[(int) box_id]++;
    }

    box_count_cummul[0] = box_count[0];
    for(unsigned n = 1; n < glob_set->natoms; n++)
    {
        box_count_cummul[n] = box_count[n] + box_count_cummul[n-1];
    }
}



void atom_box_sort()
{
    int mybox;
    box_variables_init();
    box_count_reinit();
    atom_set_box_calc();

    sotl_atom_set_t * glob_set = get_global_atom_set();

    for(unsigned i = 0; i<glob_set->natoms; i++)
    {
        mybox =atom_box_index[i];
        tmp_atom_set->pos.x[box_count_cummul[mybox]-1] = glob_set->pos.x[i];
        tmp_atom_set->pos.y[box_count_cummul[mybox]-1] = glob_set->pos.y[i];
        tmp_atom_set->pos.z[box_count_cummul[mybox]-1] = glob_set->pos.z[i];
        tmp_atom_set->speed.dx[box_count_cummul[mybox]-1] = glob_set->speed.dx[i];
        tmp_atom_set->speed.dy[box_count_cummul[mybox]-1] = glob_set->speed.dy[i];
        tmp_atom_set->speed.dz[box_count_cummul[mybox]-1] = glob_set->speed.dz[i];
        box_count_cummul[mybox]--;
    }
    for (unsigned i=0; i<glob_set->natoms; i++)
    {
        glob_set->pos.x[i] = tmp_atom_set->pos.x[i];
        glob_set->pos.y[i] = tmp_atom_set->pos.y[i];
        glob_set->pos.z[i] = tmp_atom_set->pos.z[i];
        glob_set->speed.dx[i] = tmp_atom_set->speed.dx[i];
        glob_set->speed.dy[i] = tmp_atom_set->speed.dy[i];
        glob_set->speed.dz[i] = tmp_atom_set->speed.dz[i];
    }/*
    memcpy((void *) glob_set->pos.x, tmp_atom_set->pos.x, 3*glob_set->natoms);
    memcpy((void *) glob_set->speed.dx, tmp_atom_set->speed.dx, 3*glob_set->natoms);
*/
    box_count_cummul[0]= box_count[0];
    for(unsigned n = 1; n < glob_set->natoms; n++)
    {
        box_count_cummul[n] = box_count[n] + box_count_cummul[n-1];
    }

}

void get_neighbours(index_list * neighbours)
{
    int me;
    sotl_domain_t * glob_dom = get_global_domain();
    // Choper les voisins selon x
    for(int i = 0 ; i < 1; i++)
    {
        me = neighbours->elements[i];
        neighbours->elements[neighbours->nbElements++] = me +1;
        neighbours->elements[neighbours->nbElements++] = me -1;
    }
    // Choper les voisins selon y
    for(int i = 0 ; i < 3; i++)
    {
        me = neighbours->elements[i];
        neighbours->elements[neighbours->nbElements++] = me + glob_dom->boxes[0];
        neighbours->elements[neighbours->nbElements++] = me - glob_dom->boxes[0];
    }
    for(int i = 0 ; i < 9; i++)
    {
        me = neighbours->elements[i];
        neighbours->elements[neighbours->nbElements++] = me + glob_dom->boxes[0]*glob_dom->boxes[1] ;
        neighbours->elements[neighbours->nbElements++] = me - glob_dom->boxes[0]*glob_dom->boxes[1];
    }

    /*kornflex
        for (int i =-1 ; i< 2; i++)
            for (int j =-1 ; j< 2; j++)
                for (int k =-1 ; k< 2; k++)
                {
                    if(!(i==j && j==k && k ==0))
                    {
                        neighbours->elements[neighbours->nbElements] = me + i + glob_dom->boxes[0]*(j + glob_dom->boxes[1]*k);
                        neighbours->nbElements ++;
                    }
                }
    */
    assert(neighbours-> nbElements==27);
}

// Update positions of atoms by adding (dx, dy, dz)
//
static void omp_move (sotl_device_t *dev)
{
    struct timeval tv1,tv2;
    sotl_atom_set_t *set = &dev->atom_set;

    gettimeofday(&tv1,NULL);
    #pragma omp parallel for
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
static void omp_gravity (sotl_device_t *dev)
{
    sotl_atom_set_t *set = &dev->atom_set;
    const calc_t g = 0.005;
    #pragma omp parallel for
    for (unsigned n = 0; n < set->natoms; n++)
    {
        set->speed.dz[n] -= g;
    }

}

static void omp_bounce (sotl_device_t *dev)
{
    struct timeval tv1,tv2;

    sotl_atom_set_t *set = &dev->atom_set;
    float xmin = dev->domain.min_ext[0];
    float ymin = dev->domain.min_ext[1];
    float zmin = dev->domain.min_ext[2];
    float xmax = dev->domain.max_ext[0];
    float ymax = dev->domain.max_ext[1];
    float zmax = dev->domain.max_ext[2];

    gettimeofday(&tv1,NULL);
    #pragma omp parallel for
    for (unsigned n = 0; n < set->natoms; n++)
    {
        if(set->pos.x[n] < xmin || set->pos.x[n] > xmax)
        {
            set->speed.dx[n] *= -1.0;
        }
        if(set->pos.y[n] < ymin || set->pos.y[n] > ymax)
        {
            set->speed.dy[n] *= -1.0;
        }
        if(set->pos.z[n] < zmin || set->pos.z[n] > zmax)
        {
            set->speed.dz[n] *= -1.0;
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

/*  floatcmp(float* f1, float* f2){ */
/*   return (*f1 > *f2); */
/* } */

static void omp_force (sotl_device_t *dev)
{
    struct timeval tv1,tv2;

    sotl_atom_set_t *set = &dev->atom_set;
    sotl_domain_t * glob_dom = get_global_domain();

    gettimeofday(&tv1,NULL);
    unsigned start_ind, end_ind;


    #pragma omp parallel for schedule(dynamic)
    for (unsigned current = 0; current < set->natoms; current++)
    {
        calc_t force[3] = { 0.0, 0.0, 0.0 };
        /*    index_list neighbours;
            neighbours.elements[0] = atom_box_index[current];
            //fprintf(stderr, "YOLOOO : %d\n", atom_box_index[current]);
            neighbours.current = 1;
            neighbours.nbElements=1;
            get_neighbours(&neighbours);

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
        */

        /**
          * Version Tri par boites :
          */
        /*
                for (int i = 0; i<neighbours.nbElements; i++)
                {
                    if(neighbours.elements[i] >= 0 && (unsigned)neighbours.elements[i] < glob_dom->total_boxes)
                    {
                        start_ind = box_count_cummul[neighbours.elements[i]];
                        end_ind = (glob_dom->total_boxes-1) ? set->natoms : box_count_cummul[neighbours.elements[i]+1];
                        for (unsigned other = start_ind; other<end_ind; other++)
                        {
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
                        }
                    }
                }
        */
        /** Version Tri par Z :
          */
        for (unsigned other = current+1; other < set->natoms && ((set->pos.z[other] - set->pos.z[current]) < LENNARD_CUTOFF); other++)
        {
            //Now takes the atoms near from current according to the sorted array
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
        }
        for (unsigned other = current; other > 0 && ((set->pos.z[current] - set->pos.z[other-1]) < LENNARD_CUTOFF); other--)
        {
            //Now takes the atoms near from current according to the sorted array
            if (current != other-1)
            {
                calc_t sq_dist = squared_distance (set, current, other-1);

                if (sq_dist < LENNARD_SQUARED_CUTOFF)
                {
                    calc_t intensity = lennard_jones (sq_dist);

                    force[0] += intensity * (set->pos.x[current] - set->pos.x[other-1]);
                    force[1] += intensity * (set->pos.x[set->offset + current] -
                                             set->pos.x[set->offset + other-1]);
                    force[2] += intensity * (set->pos.x[set->offset * 2 + current] -
                                             set->pos.x[set->offset * 2 + other-1]);
                }
            }
        }
        /*Fin version Tri par Z */

        set->speed.dx[current] += force[0];
        set->speed.dx[set->offset + current] += force[1];
        set->speed.dx[set->offset * 2 + current] += force[2];
    }

    gettimeofday(&tv2,NULL);
    time_force += (float)TIME_DIFF(tv1,tv2);
}


// Main simulation function
//
void omp_one_step_move (sotl_device_t *dev)
{
    // Sort positions
    sotl_atom_set_t *set = get_global_atom_set();
    for(int i = 0; i<set->natoms; i++){
        assert(set->pos.x[i] >= 0);
    }

    /* Tri par Z: */
    atom_set_sort(set);
    /* Tri par boites :*/

    // atom_box_sort();

    // Apply gravity force
    //
    if (gravity_enabled)
        omp_gravity (dev);

    // Compute interactions between atoms
    //
    if (force_enabled)
        omp_force (dev);

    // Bounce on borders
    //
    if(borders_enabled)
        omp_bounce (dev);

    // Update positions
    //
    omp_move (dev);



    n++;
    //Stops time printing after 50 iterations
    if (n<50)
        printf("\r=====time_move %.1fus - time_force %.1fus - time_bounce %.1fus=====", time_move/n, time_force/n, time_bounce/n);
    if (n==50)
        printf("\n");

    // Sort check
    //~ for (unsigned current = 0; current < set->natoms-1; current++) {
    //~ if (set->pos.z[current]> set->pos.z[current+1])
    //~ exit(0);
    //~ }

#ifdef HAVE_LIBGL
    // Update OpenGL position
    //
    if (dev->display)
        omp_update_vbo (dev);
#endif
}

void omp_init (sotl_device_t *dev)
{
#ifdef _SPHERE_MODE_
    sotl_log(ERROR, "Sequential implementation does currently not support SPHERE_MODE\n");
    exit (1);
#endif

    borders_enabled = 1;
    printf("\n");
    dev->compute = SOTL_COMPUTE_OMP; // dummy op to avoid warning

}

void omp_alloc_buffers (sotl_device_t *dev)
{
    atom_state = calloc(dev->atom_set.natoms, sizeof(int));
    printf("natoms: %d\n", dev->atom_set.natoms);
}

void omp_finalize (sotl_device_t *dev)
{
    free(atom_state);

    dev->compute = SOTL_COMPUTE_OMP; // dummy op to avoid warning
}
