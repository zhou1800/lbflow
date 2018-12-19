/*! \file file-raw.c */

#include "file-raw.h"
#include <stdio.h>

int LoadRawfile(const char *pathName, RAW *raw, int length)
{
    FILE* fp = fopen(pathName, "rb");
    VALIDATE(fp != NULL, "Can not open raw file '%s'.\n", pathName);
    VALIDATE(fread(raw, sizeof(RAW), length, fp) == length,
             "Failed to read raw file '%s'.\n", pathName);
    fclose(fp);
    return SUCCESS;
}
