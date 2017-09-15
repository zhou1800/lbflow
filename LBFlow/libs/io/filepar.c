/*! \file filepar.c
 * \author Yu Zhou (zhouyu@live.com)
 */
#include "filepar.h"
#include "config.h"
#include "cJSON.h"
#include "mpi.h"

#include <string.h>
#include <stdlib.h>

extern int ID;

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

int SyncString(char str[], int mId)
{
    int len;
    (ID == mId) && (len = strlen(str) + 1);
    /* sync string length. */
    VALIDATE(MPI_Bcast(&len, 1, MPI_INT, mId, MPI_COMM_WORLD) == MPI_SUCCESS &&
             MPI_Bcast(str, len, MPI_CHAR, mId, MPI_COMM_WORLD) == MPI_SUCCESS,
             "MPI_Bcast fail to Bcast string : '%s' on node %d.\n", str, ID);
    return NOERR;
}

int ParseParameterString(const char str[], PARAMETER* par)
{
    return NOERR;
}

int LoadParameter(const char *filePathName, PARAMETER *par)
{
    char fileContent[MAX_FILE_LEN];
    if (ID == 0 && LoadFile(filePathName, fileContent, MAX_FILE_LEN) == ERROR)
        exit(ERROR);
    return NOERR;
}
