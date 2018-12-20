/*! \file file-raw.c */

#include "file-raw.h"
#include <stdio.h>

int LoadRawfile(const char *pathName, RAW *raw, int length)
{
    int status = ERROR;

    FILE* fp = fopen(pathName, "rb");
    ASSERT(fp != NULL, "Can not open raw file '%s'.\n", pathName);
    TO_END(fread(raw, sizeof(RAW), length, fp) != length,
             "Failed to read raw file '%s'.\n", pathName);

    status = SUCCESS;
end:
    fclose(fp);
    return status;
}
