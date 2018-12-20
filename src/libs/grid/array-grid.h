/* \file array-grid.h */

#ifndef ARRAY_GRID_H
#define ARRAY_GRID_H

#include "config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief The SLICE struct is for selecting subset of a multi-dimensional array represented as
 *  an 1D array.
 */
typedef struct {
    int start;
    int sizes[MAX_DIMENSION];
    int stride[MAX_DIMENSION];
} SLICE;

/*!
 * \brief SliceToIndices calculates indices from a slice.
 * \param[in] slice - the slice.
 * \param[out] indices - array indices of the grid selected by slice.
 * \param[out] nIndices - number of indices.
 * \return SUCCESS on success, ERROR on failure.
 */
int SliceToIndices(const SLICE* slice, int* indices, int* nIndices);

/*!
 * \brief Convert an index in 1D array to it's position in grid.
 * \param[in] index - index in 1D array.
 * \param[in] size - size of grid.
 * \param[out] position - position in the grid. Must be allocated the same size of grid size.
 * \return SUCCESS on success, ERROR on failure.
 */
int IndexToPosition(int index, const int size[], int position[]);

/*!
 * \brief Convert a position in grid to it's index in 1D array.
 * \param[in] position - position in the grid.
 * \param[in] size - size of grid.
 * \param[out] index - index in 1D array.
 * \return SUCCESS on success, ERROR on failure.
 */
int PositionToIndex(const int position[], const int size[], int* index);

#ifdef __cplusplus
}
#endif

#endif /* ARRAY_GRID_H */
