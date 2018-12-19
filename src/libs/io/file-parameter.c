/*! \file file-parameter.c */

#include "file-parameter.h"
#include "cJSON.h"
#include "mpi.h"
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

extern int ID;  // Processor id

PARAMETER DEFAULT_PARAMETER = {
    {0},    // Grid
    {0},
    0,
    0,
    0,      // Control
    0,
    0,
    0,
    0,
    0,
    {0},    // Output
    {0},
    SCSP,   // Model
    SRT,
    {0},
    {0},
    0,
    0,
    0,
    0,
    {{0,0,0,0}}
};

/*!
 * \brief LoadFile load a file on disk into string.
 * \param[in] pathName - full pathe name of the file.
 * \param[out] content - file content.
 * \param[in] maxFileLength - file length must not exceed this number.
 * \return the number of characters read.
 * \return ERROR on fail to open file or incomplete file reading.
 */
int LoadFile(const char* pathName, char* content, int maxFileLength)
{
    FILE *fp = fopen(pathName, "r");
    VALIDATE(fp, "Can not open parameter file '%s'.\n", pathName);
    int fileLen = fread(content, sizeof(char), maxFileLength, fp);
    VALIDATE(fileLen < maxFileLength,
             "Incomplete file loading '%s'. "
             "File size exceeds maximum allowed size of %d bytes.\n",
             pathName, maxFileLength);
    content[fileLen] = '\0';
    fclose(fp);
    return fileLen;
}

/*!
 * \brief SyncString syncorize a string between all processors.
 * \param[in,out] str - the string to be syncorized.
 * \param[in] mId - master(source) processor id.
 * \return ERROR / SUCCESS
 */
void SyncString(char str[], int mId)
{
    int len;
    if (ID == mId) len = strlen(str) + 1;
    /* Sync string length. */
    MPI_Bcast(&len, 1, MPI_INT, mId, MPI_COMM_WORLD);
    MPI_Bcast(str, len, MPI_CHAR, mId, MPI_COMM_WORLD);
}

/*!
 * \brief ParseParameter* parse paramter in the '*' section.
 * \param[in] grid - the '*' section in JSON format.
 * \param[out] par - parse to the parameter file.
 * \return ERROR / SUCCESS
 */
int ParseGridParameter(const cJSON* grid, PARAMETER* par)
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
        VALIDATE(par->size[i] > 0, "Parameter error (gridSize).\n");
    }

    cJSON *nDirection = cJSON_GetObjectItemCaseSensitive(grid, "nDirection");
    VALIDATE(cJSON_IsNumber(nDirection), "Parameter error (nDirection).\n");
    par->nDirection = nDirection->valuedouble;

    return SUCCESS;
}

/*!
 * \see ParseGridParameter()
 */
int ParseControlParameter(const cJSON* control, PARAMETER* par)
{
    VALIDATE(control, "Parameter error (control)\n");

    cJSON *startStep = cJSON_GetObjectItemCaseSensitive(control, "startStep");
    VALIDATE(cJSON_IsNumber(startStep), "Parameter error (startStep).\n");
    par->startStep = startStep->valuedouble;

    cJSON *endStep = cJSON_GetObjectItemCaseSensitive(control, "endStep");
    VALIDATE(cJSON_IsNumber(endStep), "Parameter error (endStep).\n");
    par->endStep = endStep->valuedouble;

    cJSON *resultIntvl = cJSON_GetObjectItemCaseSensitive(control, "resultIntvl");
    VALIDATE(cJSON_IsNumber(resultIntvl), "Parameter error (resultIntvl).\n");
    par->resultIntvl = resultIntvl->valuedouble;

    cJSON *restartIntvl = cJSON_GetObjectItemCaseSensitive(control, "restartIntvl");
    VALIDATE(cJSON_IsNumber(restartIntvl), "Parameter error (restartIntvl).\n");
    par->restartIntvl = restartIntvl->valuedouble;

    cJSON *reportIntvl = cJSON_GetObjectItemCaseSensitive(control, "reportIntvl");
    VALIDATE(cJSON_IsNumber(reportIntvl), "Parameter error (reportIntvl).\n");
    par->reportIntvl = reportIntvl->valuedouble;

    cJSON *maxError = cJSON_GetObjectItemCaseSensitive(control, "maxError");
    VALIDATE(cJSON_IsNumber(maxError), "Parameter error (maxError).\n");
    par->maxError = maxError->valuedouble;

    return SUCCESS;
}

/*!
 * \see ParseGridParameter()
 */
int ParseOutputParameter(const cJSON* output, PARAMETER* par)
{
    VALIDATE(output, "Parameter error (output)\n");

    cJSON *outPath = cJSON_GetObjectItemCaseSensitive(output, "outPath");
    VALIDATE(cJSON_IsString(outPath), "Parameter error (outPath).\n");
    strcpy(par->outPath, outPath->valuestring);

    cJSON* reportROI = cJSON_GetObjectItemCaseSensitive(output, "reportROI");
    VALIDATE(cJSON_IsArray(reportROI), "Parameter error (reportROI).\n");
    int nReportROI = cJSON_GetArraySize(reportROI) < ROI_SLICE?
                     cJSON_GetArraySize(reportROI) : ROI_SLICE;
    for (int i = 0; i < nReportROI; ++i) {
        cJSON* x = cJSON_GetArrayItem(reportROI, i);
        VALIDATE(cJSON_IsNumber(x), "Parameter error (reportROI).\n");
        par->reportROI[i] = x->valuedouble;
    }
    for (int i = nReportROI; i < ROI_SLICE; ++i) {
        par->reportROI[i] = -1;
    }

    return SUCCESS;
}

/*!
 * \see ParseGridParameter()
 */
int ParseModelParameter(const cJSON* model, PARAMETER* par)
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

    cJSON *relaxModel = cJSON_GetObjectItemCaseSensitive(model, "relaxModel");
    VALIDATE(cJSON_IsString(relaxModel), "Parameter error (relaxModel).\n");
    if (!strcmp("MRT", relaxModel->valuestring))
        par->relaxModel = MRT;
    else if (!strcmp("SRT", relaxModel->valuestring))
        par->relaxModel = SRT;
    else
        VALIDATE(0, "Parameter error (relaxModel).\n");

    int nPhase = par->fluidModel == SCSP ? SINGLE_PHASE : BIN_PHASE;
    cJSON* density = cJSON_GetObjectItemCaseSensitive(model, "density");
    VALIDATE(cJSON_IsArray(density), "Parameter error (density).\n");
    VALIDATE(cJSON_GetArraySize(density) >= nPhase, "Parameter error (not enough density).\n");
    for (int i = 0; i < nPhase; ++i) {
        cJSON* x = cJSON_GetArrayItem(density, i);
        VALIDATE(cJSON_IsNumber(x), "Parameter error (density).\n");
        par->density[i] = x->valuedouble;
    }
    cJSON* relaxTime = cJSON_GetObjectItemCaseSensitive(model, "relaxTime");
    VALIDATE(cJSON_IsArray(relaxTime), "Parameter error (relaxTime).\n");
    VALIDATE(cJSON_GetArraySize(relaxTime) >= nPhase, "Parameter error (not enough relaxTime).\n");
    for (int i = 0; i < nPhase; ++i) {
        cJSON* x = cJSON_GetArrayItem(relaxTime, i);
        VALIDATE(cJSON_IsNumber(x), "Parameter error (relaxTime).\n");
        par->relaxTime[i] = x->valuedouble;
    }

    cJSON *latticeSpeed = cJSON_GetObjectItemCaseSensitive(model, "latticeSpeed");
    VALIDATE(cJSON_IsNumber(latticeSpeed), "Parameter error (latticeSpeed).\n");
    par->latticeSpeed = latticeSpeed->valuedouble;

    if (par->fluidModel == COLOR){
        cJSON *alpha0 = cJSON_GetObjectItemCaseSensitive(model, "alpha0");
        VALIDATE(cJSON_IsNumber(alpha0), "Parameter error (alpha0).\n");
        par->alpha0 = alpha0->valuedouble;

        cJSON *ift = cJSON_GetObjectItemCaseSensitive(model, "ift");
        VALIDATE(cJSON_IsNumber(ift), "Parameter error (ift).\n");
        par->ift = ift->valuedouble;

        cJSON *thickness = cJSON_GetObjectItemCaseSensitive(model, "thickness");
        VALIDATE(cJSON_IsNumber(thickness), "Parameter error (thickness).\n");
        par->thickness = thickness->valuedouble;
    }

    return SUCCESS;
}

/*!
 * \see ParseGridParameter()
 */
int ParseGrayscaleParameter(const cJSON* grayscale, PARAMETER* par)
{
    VALIDATE(grayscale, "Parameter error (grayscale)\n");

    for (int i = 0; i < cJSON_GetArraySize(grayscale); ++i){
        cJSON* item = cJSON_GetArrayItem(grayscale, i);
        // VALIDATE(item, "Parameter error (grayscale)\n");

        cJSON *index= cJSON_GetObjectItemCaseSensitive(item, "index");
        VALIDATE(cJSON_IsNumber(index), "Parameter error (index).\n");
        int indexValue = index->valuedouble;
        VALIDATE(indexValue < MAX_GRAYSCALE && indexValue >= 0,
                 "Parameter error (index %d out of range).\n", indexValue);

        cJSON *type = cJSON_GetObjectItemCaseSensitive(item, "type");
        VALIDATE(cJSON_IsString(type), "Parameter error (type of index %d).\n", indexValue);
        if (!strcmp("VD", type->valuestring)) {
            par->grayscale[indexValue].type = VD;
        }
        else if (!strcmp("GV", type->valuestring)){
            par->grayscale[indexValue].type = GV;
        }
        else if (!strcmp("PB", type->valuestring))
            par->grayscale[indexValue].type = PB;
        else if (!strcmp("VB", type->valuestring))
            par->grayscale[indexValue].type = VB;
        else if (!strcmp("SD", type->valuestring))
            par->grayscale[indexValue].type = SD;
        else
            VALIDATE(0, "Parameter error (type).\n");

        cJSON *value = cJSON_GetObjectItemCaseSensitive(item, "value");
        cJSON *direction = cJSON_GetObjectItemCaseSensitive(item, "direction");
        cJSON *phase = cJSON_GetObjectItemCaseSensitive(item, "phase");
        switch (par->grayscale[indexValue].type) {
        case VD:
            VALIDATE(cJSON_IsNumber(phase), "Parameter error (phase of index %d).\n", indexValue);
            par->grayscale[indexValue].phase = phase->valuedouble;
            break;
        case GV:
        case PB:
        case VB:
            VALIDATE(cJSON_IsNumber(value), "Parameter error (value of index %d).\n", indexValue);
            par->grayscale[indexValue].value = value->valuedouble;

            VALIDATE(cJSON_IsNumber(direction),
                     "Parameter error (direction of index %d).\n", indexValue);
            par->grayscale[indexValue].direction = direction->valuedouble;

            VALIDATE(cJSON_IsNumber(phase), "Parameter error (phase of index %d).\n", indexValue);
            par->grayscale[indexValue].phase = phase->valuedouble;
            break;
        default:    //case SD:
            VALIDATE(cJSON_IsNumber(value), "Parameter error (value of index %d).\n", indexValue);
            par->grayscale[indexValue].value = value->valuedouble;
            break;
        }
    }

    return SUCCESS;
}

/*!
 * \brief ParseParameter parse all the parameters in given string.
 * \param[in] str - parameter source string in JSON format.
 * \param[out] par - parse to the parameter file.
 * \return ERROR / SUCCESS
 */
int ParseParameter(const char str[], PARAMETER* par)
{
    int status = ERROR;

    cJSON* root = cJSON_Parse(str);
    if (!root){
        fprintf(stderr, "Parameter format error\n");
        goto end;
    }

    cJSON* grid = cJSON_GetObjectItemCaseSensitive(root, "grid");
    if (ParseGridParameter(grid, par) == ERROR)
        goto end;

    cJSON* control = cJSON_GetObjectItemCaseSensitive(root, "control");
    if (ParseControlParameter(control, par) == ERROR)
        goto end;

    cJSON* output = cJSON_GetObjectItemCaseSensitive(root, "output");
    if (ParseOutputParameter(output, par) == ERROR)
        goto end;

    cJSON* model = cJSON_GetObjectItemCaseSensitive(root, "model");
    if (ParseModelParameter(model, par) == ERROR)
        goto end;

    cJSON* grayscale = cJSON_GetObjectItemCaseSensitive(root, "grayscale");
    if (ParseGrayscaleParameter(grayscale, par) == ERROR)
        goto end;

    status = SUCCESS;

end:
    cJSON_Delete(root);
    return status;
}

int LoadParameter(const char *pathName, PARAMETER *par)
{
    char fileContent[MAX_FILE_LEN];
    if (ID == 0 && LoadFile(pathName, fileContent, MAX_FILE_LEN) == ERROR) return ERROR;
    SyncString(fileContent, 0);
    return ParseParameter(fileContent, par);
}
