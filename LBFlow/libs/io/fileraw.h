/*! \file fileraw.h
 * \author Yu Zhou (zhouyu@live.com)
 */

#ifndef FILERAW_H
#define FILERAW_H

#include "config.h"

#ifdef __cplusplus
extern "C" {
#endif

int LoadRawfile(const char* fname, RAW *raw, int size);

#ifdef __cplusplus
}
#endif

#endif // FILERAW_H
