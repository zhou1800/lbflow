/*! \file file-raw.c */

#include "file-raw.h"
#include <stdio.h>

int LoadRawfile(const char *pathName, RAW *raw, int size)
{
    FILE* fp = fopen(pathName, "rb");
    VALIDATE(fp != NULL, "Can not open raw file '%s'.\n", pathName);
    VALIDATE(fread(raw, sizeof(RAW), size, fp) == size,
             "Failed to read raw file '%s'.\n", pathName);
    fclose(fp);
    return SUCCESS;
}
