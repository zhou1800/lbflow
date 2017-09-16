/*! \file fileraw.c
 * \author Yu Zhou
 */
#include "fileraw.h"
#include <stdio.h>

int LoadRawfile(const char *fname, RAW *raw, int size)
{
    FILE* fp = fopen(fname, "rb");
    VALIDATE(fp != NULL, "Can not open raw file '%s'.\n", fname);
    // TODO: paralel reading.
    VALIDATE(fread(raw, sizeof(RAW), size, fp) == size,
                  "Failed to read raw file '%s'.\n", fname);
    fclose(fp);
    return NOERR;
}
