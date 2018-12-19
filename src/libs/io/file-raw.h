/*! \file file-raw.h */

#ifndef FILE_RAW_H
#define FILE_RAW_H

#include "config.h"

#ifdef __cplusplus
extern "C" {
#endif

int LoadRawfile(const char* pathName, RAW *raw, int size);

#ifdef __cplusplus
}
#endif

#endif /* FILE_RAW_H */
