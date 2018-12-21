/*! \file file-parameter.c */

#include "file-parameter.h"
#include "cJSON.h"
#include "mpi.h"
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

extern int ID;  /* Processor id */

PARAMETER DEFAULT_PARAMETER = {
    {0},    /* Grid */
    {0},
    0,
    0,
    0,      /* Control */
    0,
    0,
    0,
    0,
    0,
    {0},    /* Output */
    {0},
    SCSP,   /* Model */
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
 * \brief LoadFile load from a file on disk into string.
 * \param[in] pathName - full pathe name of the file.
 * \param[out] content - file content.
 * \return ERROR on fail to open file or incomplete file reading.
 */
int LoadFile(const char* pathName, char** content)
{
    FILE *fp = fopen(pathName, "r");
    ASSERT(fp, "[LoadFile] Can not open parameter file '%s'.\n", pathName);
    fseek(fp, 0L, SEEK_END);
    size_t fileLength = ftell(fp);
    rewind(fp);
    *content = (char*) malloc((fileLength + 1) * sizeof(char));
    fileLength = fread(*content, sizeof(char), fileLength, fp);
    (*content)[fileLength] = '\0';
    fclose(fp);
    return SUCCESS;
}

/*!
 * \brief SyncString syncorize a string between all processors.
 * \param[in,out] str - the string to be syncorized.
 * \param[in] mId - master(source) processor id.
 * \return ERROR / SUCCESS
 */
void SyncString(char** str, int mId)
{
    size_t len = 0;
    if (ID == mId)
        len = strlen(*str) + 1;   /* including the ending '\0' */
    /* Sync string length. */
    MPI_Bcast(&len, 1, MPI_INT, mId, MPI_COMM_WORLD);
    if (ID != mId)
        *str = (char*) malloc(len * sizeof(char));
    MPI_Bcast(*str, len, MPI_CHAR, mId, MPI_COMM_WORLD);
}

/*!
 * \brief ParseParameter* parse paramter in the '*' section.
 * \param[in] grid - the '*' section in JSON format.
 * \param[out] par - parse to the parameter file.
 * \return ERROR / SUCCESS
 */
int ParseGridParameter(const cJSON* grid, PARAMETER* par)
{
    ASSERT(grid, "[ParseGridParameter] Parameter error (grid).\n");

    cJSON *gridFile = cJSON_GetObjectItemCaseSensitive(grid, "gridFile");
    ASSERT(cJSON_IsString(gridFile), "[ParseGridParameter] Parameter error (gridFile).\n");
    strcpy(par->gridFile, gridFile->valuestring);

    cJSON* gridSize = cJSON_GetObjectItemCaseSensitive(grid, "gridSize");
    ASSERT(cJSON_IsArray(gridSize), "[ParseGridParameter] Parameter error (gridSize).\n");
    par->nDimension = cJSON_GetArraySize(gridSize) < MAX_DIMENSION ?
                      cJSON_GetArraySize(gridSize) : MAX_DIMENSION;
    for (size_t i = 0; i < par->nDimension; ++i){
        cJSON* size = cJSON_GetArrayItem(gridSize, i);
        ASSERT(cJSON_IsNumber(size), "[ParseGridParameter] Parameter error (gridSize).\n");
        par->size[i] = size->valuedouble;
        ASSERT(par->size[i] > 0, "[ParseGridParameter] Parameter error (gridSize).\n");
    }

    cJSON *nDirection = cJSON_GetObjectItemCaseSensitive(grid, "nDirection");
    ASSERT(cJSON_IsNumber(nDirection), "[ParseGridParameter] Parameter error (nDirection).\n");
    par->nDirection = nDirection->valuedouble;

    return SUCCESS;
}

/*!
 * \see ParseGridParameter()
 */
int ParseControlParameter(const cJSON* control, PARAMETER* par)
{
    ASSERT(control, "[ParseControlParameter] Parameter error (control)\n");

    cJSON *startStep = cJSON_GetObjectItemCaseSensitive(control, "startStep");
    ASSERT(cJSON_IsNumber(startStep), "[ParseControlParameter] Parameter error (startStep).\n");
    par->startStep = startStep->valuedouble;

    cJSON *endStep = cJSON_GetObjectItemCaseSensitive(control, "endStep");
    ASSERT(cJSON_IsNumber(endStep), "[ParseControlParameter] Parameter error (endStep).\n");
    par->endStep = endStep->valuedouble;

    cJSON *resultIntvl = cJSON_GetObjectItemCaseSensitive(control, "resultIntvl");
    ASSERT(cJSON_IsNumber(resultIntvl), "[ParseControlParameter] Parameter error (resultIntvl).\n");
    par->resultIntvl = resultIntvl->valuedouble;

    cJSON *restartIntvl = cJSON_GetObjectItemCaseSensitive(control, "restartIntvl");
    ASSERT(cJSON_IsNumber(restartIntvl),
           "[ParseControlParameter] Parameter error (restartIntvl).\n");
    par->restartIntvl = restartIntvl->valuedouble;

    cJSON *reportIntvl = cJSON_GetObjectItemCaseSensitive(control, "reportIntvl");
    ASSERT(cJSON_IsNumber(reportIntvl), "[ParseControlParameter] Parameter error (reportIntvl).\n");
    par->reportIntvl = reportIntvl->valuedouble;

    cJSON *maxError = cJSON_GetObjectItemCaseSensitive(control, "maxError");
    ASSERT(cJSON_IsNumber(maxError), "[ParseControlParameter] Parameter error (maxError).\n");
    par->maxError = maxError->valuedouble;

    return SUCCESS;
}

/*!
 * \see ParseGridParameter()
 */
int ParseOutputParameter(const cJSON* output, PARAMETER* par)
{
    ASSERT(output, "[ParseOutputParameter] Parameter error (output)\n");

    cJSON *outPath = cJSON_GetObjectItemCaseSensitive(output, "outPath");
    ASSERT(cJSON_IsString(outPath), "[ParseControlParameter] Parameter error (outPath).\n");
    strcpy(par->outPath, outPath->valuestring);

    cJSON* reportROI = cJSON_GetObjectItemCaseSensitive(output, "reportROI");
    ASSERT(cJSON_IsArray(reportROI), "[ParseControlParameter] Parameter error (reportROI).\n");
    int nReportROI = cJSON_GetArraySize(reportROI) < ROI_SLICE?
                     cJSON_GetArraySize(reportROI) : ROI_SLICE;
    for (int i = 0; i < nReportROI; ++i) {
        cJSON* x = cJSON_GetArrayItem(reportROI, i);
        ASSERT(cJSON_IsNumber(x), "[ParseControlParameter] Parameter error (reportROI).\n");
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
    ASSERT(model, "[ParseModelParameter] Parameter error (model)\n");

    cJSON *fluidModel = cJSON_GetObjectItemCaseSensitive(model, "fluidModel");
    ASSERT(cJSON_IsString(fluidModel), "[ParseModelParameter] Parameter error (fluidModel).\n");
    if (!strcmp("SCSP", fluidModel->valuestring))
        par->fluidModel = SCSP;
    else if (!strcmp("COLOR", fluidModel->valuestring))
        par->fluidModel = COLOR;
    else
        ASSERT(0, "[ParseModelParameter] Parameter error (fluidModel).\n");

    cJSON *relaxModel = cJSON_GetObjectItemCaseSensitive(model, "relaxModel");
    ASSERT(cJSON_IsString(relaxModel), "[ParseModelParameter] Parameter error (relaxModel).\n");
    if (!strcmp("MRT", relaxModel->valuestring))
        par->relaxModel = MRT;
    else if (!strcmp("SRT", relaxModel->valuestring))
        par->relaxModel = SRT;
    else
        ASSERT(0, "[ParseModelParameter] Parameter error (relaxModel).\n");

    int nPhase = par->fluidModel == SCSP ? SINGLE_PHASE : BIN_PHASE;
    cJSON* density = cJSON_GetObjectItemCaseSensitive(model, "density");
    ASSERT(cJSON_IsArray(density), "[ParseModelParameter] Parameter error (density).\n");
    ASSERT(cJSON_GetArraySize(density) >= nPhase,
           "[ParseModelParameter] Parameter error (not enough density).\n");
    for (int i = 0; i < nPhase; ++i) {
        cJSON* x = cJSON_GetArrayItem(density, i);
        ASSERT(cJSON_IsNumber(x), "[ParseModelParameter] Parameter error (density).\n");
        par->density[i] = x->valuedouble;
    }
    cJSON* relaxTime = cJSON_GetObjectItemCaseSensitive(model, "relaxTime");
    ASSERT(cJSON_IsArray(relaxTime), "[ParseModelParameter] Parameter error (relaxTime).\n");
    ASSERT(cJSON_GetArraySize(relaxTime) >= nPhase,
           "[ParseModelParameter] Parameter error (not enough relaxTime).\n");
    for (int i = 0; i < nPhase; ++i) {
        cJSON* x = cJSON_GetArrayItem(relaxTime, i);
        ASSERT(cJSON_IsNumber(x), "[ParseModelParameter] Parameter error (relaxTime).\n");
        par->relaxTime[i] = x->valuedouble;
    }

    cJSON *latticeSpeed = cJSON_GetObjectItemCaseSensitive(model, "latticeSpeed");
    ASSERT(cJSON_IsNumber(latticeSpeed), "[ParseModelParameter] Parameter error (latticeSpeed).\n");
    par->latticeSpeed = latticeSpeed->valuedouble;

    if (par->fluidModel == COLOR){
        cJSON *alpha0 = cJSON_GetObjectItemCaseSensitive(model, "alpha0");
        ASSERT(cJSON_IsNumber(alpha0), "[ParseModelParameter] Parameter error (alpha0).\n");
        par->alpha0 = alpha0->valuedouble;

        cJSON *ift = cJSON_GetObjectItemCaseSensitive(model, "ift");
        ASSERT(cJSON_IsNumber(ift), "[ParseModelParameter] Parameter error (ift).\n");
        par->ift = ift->valuedouble;

        cJSON *thickness = cJSON_GetObjectItemCaseSensitive(model, "thickness");
        ASSERT(cJSON_IsNumber(thickness), "[ParseModelParameter] Parameter error (thickness).\n");
        par->thickness = thickness->valuedouble;
    }

    return SUCCESS;
}

/*!
 * \see ParseGridParameter()
 */
int ParseGrayscaleParameter(const cJSON* grayscale, PARAMETER* par)
{
    ASSERT(grayscale, "[ParseGrayscaleParameter] Parameter error (grayscale)\n");

    for (int i = 0; i < cJSON_GetArraySize(grayscale); ++i){
        cJSON* item = cJSON_GetArrayItem(grayscale, i);

        cJSON *index= cJSON_GetObjectItemCaseSensitive(item, "index");
        ASSERT(cJSON_IsNumber(index), "[ParseGrayscaleParameter] Parameter error (index).\n");
        size_t indexValue = index->valuedouble;
        ASSERT(indexValue < MAX_GRAYSCALE,
               "[ParseGrayscaleParameter] Parameter error (index %ld out of range).\n", indexValue);

        cJSON *type = cJSON_GetObjectItemCaseSensitive(item, "type");
        ASSERT(cJSON_IsString(type),
               "[ParseGrayscaleParameter] Parameter error (type of index %ld).\n", indexValue);
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
            ASSERT(0, "[ParseGrayscaleParameter] Parameter error (type).\n");

        cJSON *value = cJSON_GetObjectItemCaseSensitive(item, "value");
        cJSON *direction = cJSON_GetObjectItemCaseSensitive(item, "direction");
        cJSON *phase = cJSON_GetObjectItemCaseSensitive(item, "phase");
        switch (par->grayscale[indexValue].type) {
        case VD:
            ASSERT(cJSON_IsNumber(phase),
                   "[ParseGrayscaleParameter] Parameter error (phase of index %ld).\n", indexValue);
            par->grayscale[indexValue].phase = phase->valuedouble;
            break;
        case GV:
        case PB:
        case VB:
            ASSERT(cJSON_IsNumber(value),
                   "[ParseGrayscaleParameter] Parameter error (value of index %ld).\n", indexValue);
            par->grayscale[indexValue].value = value->valuedouble;

            ASSERT(cJSON_IsNumber(direction),
                   "[ParseGrayscaleParameter] Parameter error (direction of index %ld).\n",
                   indexValue);
            par->grayscale[indexValue].direction = direction->valuedouble;

            ASSERT(cJSON_IsNumber(phase),
                   "[ParseGrayscaleParameter] Parameter error (phase of index %ld).\n", indexValue);
            par->grayscale[indexValue].phase = phase->valuedouble;
            break;
        default:    /* case SD: */
            ASSERT(cJSON_IsNumber(value),
                   "[ParseGrayscaleParameter] Parameter error (value of index %ld).\n", indexValue);
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
int ParseParameter(const char* str, PARAMETER* par)
{
    int status = ERROR;

    cJSON* root = cJSON_Parse(str);
    TO_END(!root, "[ParseParameter] Parameter format error\n");
    cJSON* grid = cJSON_GetObjectItemCaseSensitive(root, "grid");
    TO_END(ParseGridParameter(grid, par) == ERROR, "[ParseParameter]\n");
    cJSON* control = cJSON_GetObjectItemCaseSensitive(root, "control");
    TO_END(ParseControlParameter(control, par) == ERROR, "[ParseParameter]\n");
    cJSON* output = cJSON_GetObjectItemCaseSensitive(root, "output");
    TO_END(ParseOutputParameter(output, par) == ERROR, "[ParseParameter]\n")
    cJSON* model = cJSON_GetObjectItemCaseSensitive(root, "model");
    TO_END(ParseModelParameter(model, par) == ERROR, "[ParseParameter]\n")
    cJSON* grayscale = cJSON_GetObjectItemCaseSensitive(root, "grayscale");
    TO_END(ParseGrayscaleParameter(grayscale, par) == ERROR, "[ParseParameter]\n")

    status = SUCCESS;
end:
    cJSON_Delete(root);
    return status;
}

int LoadParameter(const char *pathName, PARAMETER *par)
{
    int status = ERROR;

    char* fileContent = NULL;
    int loadSuccess = ID == 0 && LoadFile(pathName, &fileContent) != ERROR;
    MPI_Bcast(&loadSuccess, 1, MPI_INT, 0, MPI_COMM_WORLD);
    TO_END(!loadSuccess, "[LoadParameter]\n");

    SyncString(&fileContent, 0);
    TO_END(ParseParameter(fileContent, par) == ERROR, "[LoadParameter]\n");

    status = SUCCESS;
end:
    free(fileContent);
    return status;
}
