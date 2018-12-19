/* \file global-grid.c */

#include "global-grid.h"
#include "config.h"
#include <mpi.h>
#include <math.h>

extern int ID;          /* Processor id */
extern int MPI_SIZE;    /* Communication size. */

/*!
 * \brief The SUBDIVISION struct defines the description of sub grid.
 */
typedef struct {
    int comm_index;         /*!< Index of communication. */
    long position[3];       /*!< location{x, y, z} of the sub grid in the global grid. */
    long size[3];           /*!< size{x, y, z} of the sub grid. */
} SUBDIVISION;

/*!
 * \brief seperate communications into dimensions, follow the image shape
 * \param[in] comm_size - number of communiations
 * \param[in] image_size_x, image_size_y, image_size_z - shape of image
 * \param[out] comm_x, comm_y, comm_y - shape of communications matrix
 */
int SeparateCommunication(
    int comm_size,
    long image_size_x, long image_size_y, long image_size_z,
    int *comm_x, int *comm_y, int *comm_z) {
    // divide processors into several dimension
    int i,j;
    long dim_size[3], swap_size;
    int dim_idx[3] = {0,1,2}, swap_idx;
    int target, dim_tgt[3], limit;
    // order dimensions, less first
    dim_size[0] = image_size_x;
    dim_size[1] = image_size_y;
    dim_size[2] = image_size_z;
    for(i = 0; i <= 1; i++) {
        for (j = i+1; j<=2; j++) {
            if (dim_size[i] > dim_size[j]) {
                swap_size = dim_size[j];
                dim_size[j] = dim_size[i];
                dim_size[i] = swap_size;
                swap_idx = dim_idx[j];
                dim_idx[j] = dim_idx[i];
                dim_idx[i] = swap_idx;
            }
        }
    }
    target = comm_size;
    for (i = 0; i<2; i++) {
        if (dim_size[i] == 1) {
            dim_tgt[dim_idx[i]] = 1;
        } else {
            dim_tgt[dim_idx[i]] = 1;
            if (i==0) {
                limit = cbrt((double) target);
            } else {
                limit = sqrt((double) target);
            }
            for (j = limit; j>1; j--) {
                if (target % j==0) {
                    dim_tgt[dim_idx[i]] = j;
                    target = target / j;
                    break;
                }
            }
        }
    }
    dim_tgt[dim_idx[2]] = target;
    *comm_x = dim_tgt[0];
    *comm_y = dim_tgt[1];
    *comm_z = dim_tgt[2];
}

int ScatterGrid(const RAW* globalGrid, const int globalSize[], RAW* localGrid,
                SLICE* localSlices, int* localSliceIndex) {
    return SUCCESS;
}
