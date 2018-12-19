/* \file global-grid.h */

#ifndef GLOBAL_GRID_H
#define GLOBAL_GRID_H

#include "config.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int start;
    int sizes[MAX_DIMENSION];
    int stride[MAX_DIMENSION];
} SLICE;

/*!
 * \brief ScatterGrid divide global grid into sub grids, and send them to slave cores.
 * \param[in] comm_size - number of communiations
 * \param[in] image_size_x, image_size_y, image_size_z - shape of image
 * \param[out] comm_x, comm_y, comm_y - shape of communications matrix
 */
int ScatterGrid(const RAW* globalGrid, const int globalSize[], RAW* localGrid,
                SLICE* localSlices, int* localSliceIndex);

#ifdef __cplusplus
}
#endif

#endif // GLOBAL_GRID_H
