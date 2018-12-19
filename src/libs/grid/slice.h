/* \file slice.h */

#ifndef SLICE_H
#define SLICE_H

#include "config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief Convert an index in 1D array to it's location in grid.
 * \param[in] index - index in 1D array.
 * \param[in] size - size of grid.
 * \param[out] location - location in the grid. Must be allocated the same size of grid size.
 * \return SUCCESS on success, ERROR on failure.
 */
int IndexToLocation(int index, const int size[], int location[]);

/*!
 * \brief Convert a location in grid to it's index in 1D array.
 * \param[in] location - location in the grid.
 * \param[in] size - size of grid.
 * \param[out] index - index in 1D array.
 * \return SUCCESS on success, ERROR on failure.
 */
int LocationToIndex(const int location[], const int size[], int* index);

#ifdef __cplusplus
}
#endif

#endif /* SLICE_H */
