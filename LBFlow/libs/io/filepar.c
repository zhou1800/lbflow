/*! \file filepar.c
 * \author Yu Zhou (zhouyu@live.com)
 */
#include "filepar.h"
#include "config.h"
#include "cJSON.h"
#include "mpi.h"
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

extern int ID;

PARAMETER DEFAULT_PARAMETER = {
    {0},    // grid
    {0},
    0,
    0,
    0,
    0,      // control
    0,
    0,
    0,
    0,
    0,
    {0},    // output
    0,
    {0},
    SCSP,      // fluid
    SRT,
    {0},
    0,
    {0},
    0,
    0,
    0,
    {{0,0,0,0}}
};

/*!
 * \brief LoadFile load a file on disk into string.
 * \param filePathName[in] - full pathe name of the file.
 * \param content[out] - file content.
 * \param maxFileLength[in] - file length must not exceed this number.
 * \return the number of characters read.
 * \return ERROR on fail to open file or incomplete file reading.
 */
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

/*!
 * \brief SyncString syncorize a string between all processors.
 * \param str[in,out] - the string to be syncorized.
 * \param mId[in] - master(source) processor id.
 * \return ERROR / NOERR
 */
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

/*!
 * \brief ParseParameter* parse paramter in the '*' section.
 * \param grid[in] - the '*' section in JSON format.
 * \param par[out] - parse to the parameter file.
 * \return ERROR / NOERR
 */
int ParseParameterGrid(const cJSON* grid, PARAMETER* par)
{
    VALIDATE(grid, "Parameter error (grid).\n");

    cJSON *gridFile = cJSON_GetObjectItemCaseSensitive(grid, "gridFile");
    VALIDATE(cJSON_IsString(gridFile), "Parameter error (gridFile).\n");
    strcpy(par->gridFile, gridFile->valuestring);

    cJSON* gridSize = cJSON_GetObjectItemCaseSensitive(grid, "gridSize");
    VALIDATE(cJSON_IsArray(gridSize), "Parameter error (gridSize).\n");
    par->nDimension = cJSON_GetArraySize(gridSize) < MAX_DIMENSION ?
                cJSON_GetArraySize(gridSize) : MAX_DIMENSION;
    for (int i = 0; i < par->nDimension; ++i){
        cJSON* size = cJSON_GetArrayItem(gridSize, i);
        VALIDATE(cJSON_IsNumber(size), "Parameter error (gridSize).\n");
        par->size[i] = size->valuedouble;
    }

    cJSON *nDirection = cJSON_GetObjectItemCaseSensitive(grid, "nDirection");
    VALIDATE(cJSON_IsNumber(nDirection), "Parameter error (nDirection).\n");
    par->nDirection = nDirection->valuedouble;

    cJSON *latticeLen = cJSON_GetObjectItemCaseSensitive(grid, "latticeLen");
    VALIDATE(cJSON_IsNumber(latticeLen), "Parameter error (latticeLen).\n");
    par->latticeLen = latticeLen->valuedouble;

    return NOERR;
}

/*!
 * \see ParseParameterGrid()
 */
int ParseParameterControl(const cJSON* control, PARAMETER* par)
{
    VALIDATE(control, "Parameter error (control)\n");

    cJSON *startTime = cJSON_GetObjectItemCaseSensitive(control, "startTime");
    VALIDATE(cJSON_IsNumber(startTime), "Parameter error (startTime).\n");
    par->startTime = startTime->valuedouble;

    cJSON *endTime = cJSON_GetObjectItemCaseSensitive(control, "endTime");
    VALIDATE(cJSON_IsNumber(endTime), "Parameter error (endTime).\n");
    par->endTime = endTime->valuedouble;

    cJSON *writeIntvl = cJSON_GetObjectItemCaseSensitive(control, "writeIntvl");
    VALIDATE(cJSON_IsNumber(writeIntvl), "Parameter error (writeIntvl).\n");
    par->writeIntvl = writeIntvl->valuedouble;

    cJSON *rptIntvl = cJSON_GetObjectItemCaseSensitive(control, "rptIntvl");
    VALIDATE(cJSON_IsNumber(rptIntvl), "Parameter error (rptIntvl).\n");
    par->rptIntvl = rptIntvl->valuedouble;

    cJSON *rstIntvl = cJSON_GetObjectItemCaseSensitive(control, "rstIntvl");
    VALIDATE(cJSON_IsNumber(rstIntvl), "Parameter error (rstIntvl).\n");
    par->rstIntvl = rstIntvl->valuedouble;

    cJSON *maxError = cJSON_GetObjectItemCaseSensitive(control, "maxError");
    VALIDATE(cJSON_IsNumber(maxError), "Parameter error (maxError).\n");
    par->maxError = maxError->valuedouble;

    return NOERR;
}

/*!
 * \see ParseParameterGrid()
 */
int ParseParameterOutput(const cJSON* output, PARAMETER* par)
{
    VALIDATE(output, "Parameter error (output)\n");

    cJSON *outPath = cJSON_GetObjectItemCaseSensitive(output, "outPath");
    VALIDATE(cJSON_IsString(outPath), "Parameter error (outPath).\n");
    strcpy(par->outPath, outPath->valuestring);

    cJSON *rptDir = cJSON_GetObjectItemCaseSensitive(output, "rptDir");
    VALIDATE(cJSON_IsNumber(rptDir), "Parameter error (rptDir).\n");
    par->rptDir = rptDir->valuedouble;

    if (par->rptDir != 0) {
        cJSON* rptSlice = cJSON_GetObjectItemCaseSensitive(output, "rptSlice");
        VALIDATE(cJSON_IsArray(rptSlice), "Parameter error (rptSlice).\n");
        VALIDATE(cJSON_GetArraySize(rptSlice) >= NUMBER_OF_SLICE,
                 "Parameter error (not enough rptSlice).\n");
        for (int i = 0; i < NUMBER_OF_SLICE; ++i) {
            cJSON* slice = cJSON_GetArrayItem(rptSlice, i);
            VALIDATE(cJSON_IsNumber(slice), "Parameter error (rptSlice).\n");
            par->rptSlice[i] = slice->valuedouble;
        }
    }

    return NOERR;
}

/*!
 * \see ParseParameterGrid()
 */
int ParseParameterModel(cJSON* model, PARAMETER* par)
{
    VALIDATE(model, "Parameter error (model)\n");

    cJSON *fluidModel = cJSON_GetObjectItemCaseSensitive(model, "fluidModel");
    VALIDATE(cJSON_IsString(fluidModel), "Parameter error (fluidModel).\n");
    if (!strcmp("SCSP", fluidModel->valuestring))
        par->fluidModel = SCSP;
    else if (!strcmp("COLOR", fluidModel->valuestring))
        par->fluidModel = COLOR;
    else
        VALIDATE(0, "Parameter error (fluidModel).\n");

    cJSON *equilModel = cJSON_GetObjectItemCaseSensitive(model, "equilModel");
    VALIDATE(cJSON_IsString(equilModel), "Parameter error (equilModel).\n");
    if (!strcmp("MRT", equilModel->valuestring))
        par->equilModel = MRT;
    else if (!strcmp("SRT", equilModel->valuestring))
        par->equilModel = SRT;
    else
        VALIDATE(0, "Parameter error (equilModel).\n");

    cJSON* rlxTime = cJSON_GetObjectItemCaseSensitive(model, "rlxTime");
    VALIDATE(cJSON_IsArray(rlxTime), "Parameter error (rlxTime).\n");
    cJSON* density = cJSON_GetObjectItemCaseSensitive(model, "density");
    VALIDATE(cJSON_IsArray(density), "Parameter error (density).\n");
    int nPhase = 0;
    if (par->fluidModel == SCSP) {
        VALIDATE(cJSON_GetArraySize(rlxTime) >= SINGLE_PHASE,
                 "Parameter error (not enough rlxTime).\n");
        VALIDATE(cJSON_GetArraySize(density) >= SINGLE_PHASE,
                 "Parameter error (not enough density).\n");
        nPhase = SINGLE_PHASE;
    }
    else if (par->fluidModel == COLOR){
        VALIDATE(cJSON_GetArraySize(rlxTime) >= BIN_PHASE,
                 "Parameter error (not enough rlxTime).\n");
        VALIDATE(cJSON_GetArraySize(density) >= BIN_PHASE,
                 "Parameter error (not enough density).\n");
        nPhase = BIN_PHASE;
    }
    for (int i = 0; i < nPhase; ++i) {
        cJSON* slice = cJSON_GetArrayItem(rlxTime, i);
        VALIDATE(cJSON_IsNumber(slice), "Parameter error (rlxTime).\n");
        par->rlxTime[i] = slice->valuedouble;
        cJSON* dens = cJSON_GetArrayItem(density, i);
        VALIDATE(cJSON_IsNumber(dens), "Parameter error (density).\n");
        par->density[i] = dens->valuedouble;
    }

    cJSON *latticeSpd = cJSON_GetObjectItemCaseSensitive(model, "latticeSpd");
    VALIDATE(cJSON_IsNumber(latticeSpd), "Parameter error (latticeSpd).\n");
    par->latticeSpd = latticeSpd->valuedouble;

    if (par->fluidModel == COLOR){
        cJSON *alpha0 = cJSON_GetObjectItemCaseSensitive(model, "alpha0");
        VALIDATE(cJSON_IsNumber(alpha0), "Parameter error (alpha0).\n");
        par->alpha0 = alpha0->valuedouble;

        cJSON *ift = cJSON_GetObjectItemCaseSensitive(model, "ift");
        VALIDATE(cJSON_IsNumber(ift), "Parameter error (ift).\n");
        par->ift = ift->valuedouble;

        cJSON *intfThk = cJSON_GetObjectItemCaseSensitive(model, "intfThk");
        VALIDATE(cJSON_IsNumber(intfThk), "Parameter error (intfThk).\n");
        par->intfThk = intfThk->valuedouble;
    }

    return NOERR;
}

/*!
 * \see ParseParameterGrid()
 */
int ParseParameterIndex(cJSON* index, PARAMETER* par)
{
    VALIDATE(index, "Parameter error (index)\n");

    for (int i = 0; i < cJSON_GetArraySize(index); ++i){
        cJSON* item = cJSON_GetArrayItem(index, i);
        VALIDATE(item, "Parameter error (index)\n");

        cJSON *idx = cJSON_GetObjectItemCaseSensitive(item, "idx");
        VALIDATE(cJSON_IsNumber(idx), "Parameter error (idx).\n");
        int idxValue = idx->valuedouble;
        VALIDATE(idxValue < MAX_GRAY && idxValue >= 0,
                 "Parameter error (idx %d out of range).\n", idxValue);

        cJSON *type = cJSON_GetObjectItemCaseSensitive(item, "type");
        VALIDATE(cJSON_IsString(type),
                 "Parameter error (type of idx %d).\n", idxValue);
        if (!strcmp("VD", type->valuestring)) {
            par->index[idxValue].type = VD;
        }
        else if (!strcmp("GV", type->valuestring)){
            par->index[idxValue].type = GV;
        }
        else if (!strcmp("PB", type->valuestring))
            par->index[idxValue].type = PB;
        else if (!strcmp("VB", type->valuestring))
            par->index[idxValue].type = VB;
        else if (!strcmp("SD", type->valuestring))
            par->index[idxValue].type = SD;
        else
            VALIDATE(0, "Parameter error (type).\n");

        cJSON *value = cJSON_GetObjectItemCaseSensitive(item, "value");
        cJSON *dir = cJSON_GetObjectItemCaseSensitive(item, "dir");
        cJSON *phase = cJSON_GetObjectItemCaseSensitive(item, "phase");
        switch (par->index[idxValue].type) {
        case VD:
            VALIDATE(cJSON_IsNumber(phase),
                     "Parameter error (phase of idx %d).\n", idxValue);
            par->index[idxValue].phase = phase->valuedouble;
            break;
        case GV:
        case PB:
        case VB:
            VALIDATE(cJSON_IsNumber(value),
                     "Parameter error (value of idx %d).\n", idxValue);
            par->index[idxValue].value = value->valuedouble;

            VALIDATE(cJSON_IsNumber(dir),
                     "Parameter error (dir of idx %d).\n", idxValue);
            par->index[idxValue].dir = dir->valuedouble;

            VALIDATE(cJSON_IsNumber(phase),
                     "Parameter error (phase of idx %d).\n", idxValue);
            par->index[idxValue].phase = phase->valuedouble;
            break;
        case SD:
            VALIDATE(cJSON_IsNumber(value),
                     "Parameter error (value of idx %d).\n", idxValue);
            par->index[idxValue].value = value->valuedouble;
            break;
        }
    }

    return NOERR;
}

/*!
 * \brief ParseParameter parse all the parameters in given string.
 * \param str[in] - parameter source string in JSON format.
 * \param par[out] - parse to the parameter file.
 * \return ERROR / NOERR
 */
int ParseParameter(const char str[], PARAMETER* par)
{
    cJSON* root = cJSON_Parse(str);
    VALIDATE(root, "Parameter format error\n");

    cJSON* grid = cJSON_GetObjectItemCaseSensitive(root, "grid");
    if (ParseParameterGrid(grid, par) == ERROR)
        return ERROR;

    cJSON* control = cJSON_GetObjectItemCaseSensitive(root, "control");
    if (ParseParameterControl(control, par) == ERROR)
        return ERROR;

    cJSON* output = cJSON_GetObjectItemCaseSensitive(root, "output");
    if (ParseParameterOutput(output, par) == ERROR)
        return ERROR;

    cJSON* model = cJSON_GetObjectItemCaseSensitive(root, "model");
    if (ParseParameterModel(model, par) == ERROR)
        return ERROR;

    cJSON* index = cJSON_GetObjectItemCaseSensitive(root, "index");
    if (ParseParameterIndex(index, par) == ERROR)
        return ERROR;

    cJSON_Delete(root);
    return NOERR;
}

int LoadParameter(const char *filePathName, PARAMETER *par)
{
    char fileContent[MAX_FILE_LEN];
    if (ID == 0 && LoadFile(filePathName, fileContent, MAX_FILE_LEN) == ERROR)
        exit(ERROR);
    if (SyncString(fileContent, 0) == ERROR)
        exit(ERROR);
    if (ParseParameter(fileContent, par) == ERROR)
        exit(ERROR);
    return NOERR;
}
