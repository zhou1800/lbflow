/* \file global-grid.h */

#ifndef GLOBAL_GRID_H
#define GLOBAL_GRID_H

#include "array-grid.h"
#include "config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief ScatterGrid divide global grid into sub grids, and send them to slave cores.
 * \param[in] globalGrid - Raw data of global grid (required for master, NULL for slaves).
 * \param[in] globalSize - Global grid size (required for master, NULL for slaves).
 * \param[out] localGrid - Raw data for local grid.
 * \param[out] localSlices - Array of local Slices.
 * \return SUCCESS on success, ERROR on failure.
 */
int ScatterGrid(const RAW* globalGrid, const size_t globalSize[], RAW* localGrid, SLICE* localSlices);

#ifdef __cplusplus
}
#endif

#endif /* GLOBAL_GRID_H */
