/* \file global-grid.h */

#ifndef GLOBAL_GRID_H
#define GLOBAL_GRID_H

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief get the total processors of mpi.
 */
int GetCommunicationSize();

int SeparateCommunication(
    int comm_size,
    long image_size_x, long image_size_y, long image_size_z,
    int *comm_x, int *comm_y, int *comm_z);

#ifdef __cplusplus
}
#endif

#endif // GLOBAL_GRID_H
