/* \file global-grid.h */

#ifndef GLOBAL_GRID_H
#define GLOBAL_GRID_H

<<<<<<< HEAD
=======
#include "config.h"

>>>>>>> c0ac77790f9a604db6032bc86009a58fd27bee3e
#ifdef __cplusplus
extern "C" {
#endif

<<<<<<< HEAD
/*!
 * \brief The SUBDIVISION struct defines the description of sub grid.
 */
typedef struct {
    int comm_index;         /*!< Index of communication. */
    long position[3];       /*!< location{x, y, z} of the sub grid in the global grid. */
    long size[3];           /*!< size{x, y, z} of the sub grid. */
} SUBDIVISION;

/*!
 * \brief get the total processors of mpi.
 */
int GetCommunicationSize();

/*!
 * \brief seperate communications into dimensions, follow the image shape
=======
typedef struct {
    int start;
    int sizes[MAX_DIMENSION];
    int stride[MAX_DIMENSION];
} SLICE;

/*!
 * \brief ScatterGrid divide global grid into sub grids, and send them to slave cores.
>>>>>>> c0ac77790f9a604db6032bc86009a58fd27bee3e
 * \param[in] comm_size - number of communiations
 * \param[in] image_size_x, image_size_y, image_size_z - shape of image
 * \param[out] comm_x, comm_y, comm_y - shape of communications matrix
 */
<<<<<<< HEAD
int SeparateCommunication(
    int comm_size,
    long image_size_x, long image_size_y, long image_size_z,
    int *comm_x, int *comm_y, int *comm_z);
=======
int ScatterGrid(const RAW* globalGrid, const int globalSize[], RAW* localGrid,
                SLICE* localSlices, int* localSliceIndex);
>>>>>>> c0ac77790f9a604db6032bc86009a58fd27bee3e

#ifdef __cplusplus
}
#endif

#endif // GLOBAL_GRID_H
