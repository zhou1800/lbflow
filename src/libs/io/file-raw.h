/*! \file file-raw.h */

#ifndef FILE_RAW_H
#define FILE_RAW_H

#include "config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief The LoadParameter function loads the parameter file, and then parse it.
 *
 * \param[in] pathName - parameter file path name.
 * \param[out] raw - raw data from loaded file.
 * \param[out] length - raw data length (should be smaller or equal to file length).
 * \return SUCCESS on success, ERROR on failure.
 */
int LoadRawfile(const char* pathName, RAW *raw, int length);

#ifdef __cplusplus
}
#endif

#endif /* FILE_RAW_H */
