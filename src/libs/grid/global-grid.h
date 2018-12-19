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

/**!
 * \brief seperate communications into dimensions, follow the image shape
 * \param[in] comm_size - number of communiations
 * \param[in] image_size_x, image_size_y, image_size_z - shape of image
 * \param[out] comm_x, comm_y, comm_y - shape of communications matrix
 */
int SeparateCommunication(
    int comm_size,
    long image_size_x, long image_size_y, long image_size_z,
    int *comm_x, int *comm_y, int *comm_z);

#ifdef __cplusplus
}
#endif

#endif // GLOBAL_GRID_H
