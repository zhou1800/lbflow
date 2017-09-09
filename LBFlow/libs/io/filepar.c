/*! \file filepar.c
 * \author Yu Zhou
 */
#include "filepar.h"
#include "config.h"
#include "cJSON.h"

int LoadFile(const char* filePathName, char* content, int maxFileLength)
{
    FILE *fp = fopen(filePathName, "r");
    VALIDATE(fp, "Can not open parameter file '%s'.\n", filePathName);
    int fileLen = fread(content, sizeof(char), maxFileLength, fp);
    content[fileLen] = '\0';
    VALIDATE(fileLen < maxFileLength,
             "(Possible) incomplete file loading '%s'. "
             "File size exceeds maximum allowed size of %d bytes.\n",
             filePathName, maxFileLength);
    fclose(fp);
    return fileLen;
}

int LoadParameter(const char *filePathName, PARAMETER *par)
{
    return NOERR;
}
