/* \file grid-common.h */

#ifndef GRID_COMMON_H
#define GRID_COMMON_H

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief Initialise special related global variable.
 * \param[in] dimension - Number of dimensions.
 * \return SUCCESS on success, ERROR on failure.
 */
int Init(int dimension);

#ifdef __cplusplus
}
#endif

#endif /* GRID_COMMON_H */
