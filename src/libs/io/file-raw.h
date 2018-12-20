/*! \file file-raw.h */

#ifndef FILE_RAW_H
#define FILE_RAW_H

#include "config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief The LoadRawFile function loads the raw file.
 *
 * \param[in] pathName - parameter file path name.
 * \param[out] raw - raw data from loaded file.
 * \param[in] length - raw data length.
 * \return SUCCESS on success, ERROR on failure.
 */
int LoadRawFile(const char* pathName, RAW** raw, int length);

#ifdef __cplusplus
}
#endif

#endif /* FILE_RAW_H */
