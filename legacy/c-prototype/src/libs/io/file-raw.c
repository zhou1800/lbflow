/*! \file file-raw.c */

#include "file-raw.h"
#include <stdio.h>
#include <stdlib.h>

int LoadRawFile(const char* pathName, RAW** raw, size_t length)
{
    int status = ERROR;

    FILE* fp = fopen(pathName, "rb");
    ASSERT(fp != NULL, "[LoadRawFile] Can not open raw file '%s'.\n", pathName);
    *raw = (RAW*) malloc(length * sizeof(RAW));
    TO_END(fread(*raw, sizeof(RAW), length, fp) != length,
           "[LoadRawFile] Failed to read raw file '%s'.\n", pathName);

    status = SUCCESS;
end:
    fclose(fp);
    return status;
}
