/*! \file file-parameter-test.h */

#ifndef FILE_PARAMETER_TEST_H
#define FILE_PARAMETER_TEST_H

#include "file-parameter.h"
#include "../test-cfg.h"
#include "../gtest-helper.h"
#include "gtest/gtest.h"
#include "cJSON.h"

extern "C" {
    extern int ID;
    extern PARAMETER DEFAULT_PARAMETER;
    int LoadFile(const char* filePathName, char** content);
    int SyncString(char** str, int mId);
    int ParseGridParameter(const cJSON* grid, PARAMETER* par);
    int ParseControlParameter(const cJSON* control, PARAMETER* par);
    int ParseOutputParameter(const cJSON* output, PARAMETER* par);
    int ParseModelParameter(const cJSON* model, PARAMETER* par);
    int ParseGrayscaleParameter(const cJSON* grayscale, PARAMETER* par);
    int ParseParameter(const char str[], PARAMETER* par);
}

struct PARAMETER_FILE{
    const char* parameterFile;
    const char* content;
    PARAMETER par;
};

PARAMETER_FILE Grid2D{
    SOURCE_DIR"examples/Grid2D.json",     // parameter file
    "{\n"
    "    \"grid\" : {\n"
    "        \"gridFile\":     \"Grid20x12y.raw\",\n"
    "        \"gridSize\":     [20, 12],\n"
    "        \"nDirection\":   9\n"
    "    },\n"
    "    \"control\" : {\n"
    "        \"startStep\":    0,\n"
    "        \"endStep\":      1000,\n"
    "        \"resultIntvl\":  100,\n"
    "        \"restartIntvl\": 500,\n"
    "        \"reportIntvl\":  10,\n"
    "        \"maxError\":     0.000001\n"
    "    },\n"
    "    \"output\" : {\n"
    "        \"outPath\":      \"./LBFlow-Result/Grid20x12y\",\n"
    "        \"reportROI\":    [0, 98]\n"
    "    },\n"
    "    \"model\" : {\n"
    "        \"fluidModel\":   \"SCSP\",\n"
    "        \"relaxModel\":   \"SRT\",\n"
    "        \"density\":      [1],\n"
    "        \"relaxTime\":    [1],\n"
    "        \"latticeSpeed\": 1\n"
    "    },\n"
    "    \"grayscale\" : [\n"
    "        {\n"
    "            \"index\": 0,\n"
    "            \"type\": \"VD\",\n"
    "            \"phase\": 0\n"
    "        },\n"
    "        {\n"
    "            \"index\": 50,\n"
    "            \"type\": \"VB\",\n"
    "            \"value\": 0.00001,\n"
    "            \"direction\": 1,\n"
    "            \"phase\": 0\n"
    "        },\n"
    "        {\n"
    "            \"index\": 75,\n"
    "            \"type\": \"GV\",\n"
    "            \"value\": 0.001,\n"
    "            \"direction\": 1,\n"
    "            \"phase\": 0\n"
    "        },\n"
    "        {\n"
    "            \"index\": 100,\n"
    "            \"type\": \"VD\",\n"
    "            \"phase\": 1\n"
    "        },\n"
    "        {\n"
    "            \"index\": 150,\n"
    "            \"type\": \"VB\",\n"
    "            \"value\": 0.00001,\n"
    "            \"direction\": 1,\n"
    "            \"phase\": 1\n"
    "        },\n"
    "        {\n"
    "            \"index\": 200,\n"
    "            \"type\": \"PB\",\n"
    "            \"value\": 1,\n"
    "            \"direction\": 2,\n"
    "            \"phase\": 0\n"
    "        },\n"
    "        {\n"
    "            \"index\": 250,\n"
    "            \"type\": \"SD\",\n"
    "            \"value\": 90\n"
    "        },\n"
    "        {\n"
    "            \"index\": 255,\n"
    "            \"type\": \"SD\",\n"
    "            \"value\": 160\n"
    "        }\n"
    "    ]\n"
    "}\n",
    // parameters
    {
        "Grid20x12y.raw",
        {20, 12},
        2,
        9,
        0,
        1000,
        100,
        500,
        10,
        0.000001,
        "./LBFlow-Result/Grid20x12y",
        {0, 98, -1, -1, -1, -1},
        SCSP,
        SRT,
        {1, 0},
        {1, 0},
        1,
        0,
        0,
        0,
        {{VD, 0, 0, 0}, {}, {}, {}, {}, {}, {}, {}, {}, {},         //0-9
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //10-19
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //20-29
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //30-39
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //40-49
         {VB, 0.00001, 1, 0}, {}, {}, {}, {}, {}, {}, {}, {}, {},   //50-59
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //60-69
         {}, {}, {}, {}, {}, {GV, 0.001, 1, 0}, {}, {}, {}, {},     //70-79
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //80-89
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //90-99
         {VD, 0, 0, 1}, {}, {}, {}, {}, {}, {}, {}, {}, {},         //100-109
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //110-119
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //120-129
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //130-139
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //140-149
         {VB, 0.00001, 1, 1}, {}, {}, {}, {}, {}, {}, {}, {}, {},   //150-159
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //160-169
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //170-179
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //180-189
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //190-199
         {PB, 1, 2, 0}, {}, {}, {}, {}, {}, {}, {}, {}, {},         //200-209
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //210-219
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //220-229
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //230-239
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //240-249
         {SD, 90, 0, 0}, {}, {}, {}, {}, {SD, 160, 0, 0},           //250-255
        }
    },
};

PARAMETER_FILE Channel3D{
    SOURCE_DIR"examples/Channel3D.json",  // parameter file
    "{\n"
    "    \"grid\" : {\n"
    "        \"gridFile\":     \"Channel3D-24x24y100z.raw\",\n"
    "        \"gridSize\":     [24, 24, 100],\n"
    "        \"nDirection\":   19\n"
    "    },\n"
    "    \"control\" : {\n"
    "        \"startStep\":    0,\n"
    "        \"endStep\":      1000,\n"
    "        \"resultIntvl\":  100,\n"
    "        \"restartIntvl\": 500,\n"
    "        \"reportIntvl\":  10,\n"
    "        \"maxError\":     0.000001\n"
    "    },\n"
    "    \"output\" : {\n"
    "        \"outPath\":      \"./LBFlow-Result/Channel3D-24x24y100z\",\n"
    "        \"reportROI\":    [-1, -1, -1, -1, 98, -1, 77]\n"
    "    },\n"
    "    \"model\" : {\n"
    "        \"fluidModel\":   \"COLOR\",\n"
    "        \"relaxModel\":   \"MRT\",\n"
    "        \"density\":      [1, 0.831683168],\n"
    "        \"relaxTime\":    [0.57, 1.215416667],\n"
    "        \"latticeSpeed\": 1,\n"
    "        \"alpha0\":       0.33333333333333333333,\n"
    "        \"ift\":          0.085920139,\n"
    "        \"thickness\":    0.7\n"
    "    },\n"

    "    \"grayscale\" : [\n"
    "        {\n"
    "            \"index\": 0,\n"
    "            \"type\": \"VD\",\n"
    "            \"phase\": 0\n"
    "        },\n"
    "        {\n"
    "            \"index\": 50,\n"
    "            \"type\": \"VB\",\n"
    "            \"value\": 0.00001,\n"
    "            \"direction\": 1,\n"
    "            \"phase\": 0\n"
    "        },\n"
    "        {\n"
    "            \"index\": 75,\n"
    "            \"type\": \"GV\",\n"
    "            \"value\": 0.001,\n"
    "            \"direction\": 1,\n"
    "            \"phase\": 0\n"
    "        },\n"
    "        {\n"
    "            \"index\": 100,\n"
    "            \"type\": \"VD\",\n"
    "            \"phase\": 1\n"
    "        },\n"
    "        {\n"
    "            \"index\": 150,\n"
    "            \"type\": \"VB\",\n"
    "            \"value\": 0.00001,\n"
    "            \"direction\": 1,\n"
    "            \"phase\": 1\n"
    "        },\n"
    "        {\n"
    "            \"index\": 200,\n"
    "            \"type\": \"PB\",\n"
    "            \"value\": 1,\n"
    "            \"direction\": 2,\n"
    "            \"phase\": 0\n"
    "        },\n"
    "        {\n"
    "            \"index\": 250,\n"
    "            \"type\": \"SD\",\n"
    "            \"value\": 90\n"
    "        },\n"
    "        {\n"
    "            \"index\": 255,\n"
    "            \"type\": \"SD\",\n"
    "            \"value\": 160\n"
    "        }\n"
    "    ]\n"
    "}\n",
    {
        "Channel3D-24x24y100z.raw",
        {24, 24, 100},
        3,
        19,
        0,
        1000,
        100,
        500,
        10,
        0.000001,
        "./LBFlow-Result/Channel3D-24x24y100z",
        {-1, -1, -1, -1, 98, -1},
        COLOR,
        MRT,
        {1, 0.831683168},
        {0.57, 1.215416667},
        1,
        0.33333333333333333333,
        0.085920139,
        0.7,
        {{VD, 0, 0, 0}, {}, {}, {}, {}, {}, {}, {}, {}, {},         //0-9
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //10-19
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //20-29
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //30-39
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //40-49
         {VB, 0.00001, 1, 0}, {}, {}, {}, {}, {}, {}, {}, {}, {},   //50-59
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //60-69
         {}, {}, {}, {}, {}, {GV, 0.001, 1, 0}, {}, {}, {}, {},     //70-79
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //80-89
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //90-99
         {VD, 0, 0, 1}, {}, {}, {}, {}, {}, {}, {}, {}, {},         //100-109
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //110-119
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //120-129
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //130-139
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //140-149
         {VB, 0.00001, 1, 1}, {}, {}, {}, {}, {}, {}, {}, {}, {},   //150-159
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //160-169
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //170-179
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //180-189
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //190-199
         {PB, 1, 2, 0}, {}, {}, {}, {}, {}, {}, {}, {}, {},         //200-209
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //210-219
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //220-229
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //230-239
         {}, {}, {}, {}, {}, {}, {}, {}, {}, {},                    //240-249
         {SD, 90, 0, 0}, {}, {}, {}, {}, {SD, 160, 0, 0},           //250-255
        }
    },
};

class FileParameterTest : public ::testing::TestWithParam<PARAMETER_FILE>{
protected:
    virtual void SetUp(){
        fileContent = NULL;
    }
    virtual void TearDown(){
        free(fileContent);
    }
    char* fileContent;
};

INSTANTIATE_TEST_CASE_P(2Dand3D, FileParameterTest, ::testing::Values(Grid2D, Channel3D));

void CompareParameter(const PARAMETER& par1, const PARAMETER& par2)
{
    EXPECT_STREQ(par1.gridFile, par2.gridFile);
    EXPECT_ARR_EQ(par1.size, par2.size, MAX_DIMENSION);
    EXPECT_EQ(par1.nDimension, par2.nDimension);
    EXPECT_EQ(par1.nDirection, par2.nDirection);
    EXPECT_EQ(par1.startStep, par2.startStep);
    EXPECT_EQ(par1.endStep, par2.endStep);
    EXPECT_EQ(par1.resultIntvl, par2.resultIntvl);
    EXPECT_EQ(par1.restartIntvl, par2.restartIntvl);
    EXPECT_EQ(par1.reportIntvl, par2.reportIntvl);
    EXPECT_FLOAT_EQ(par1.maxError, par2.maxError);
    EXPECT_STREQ(par1.outPath, par2.outPath);
    EXPECT_ARR_EQ(par1.reportROI, par2.reportROI, ROI_SLICE);
    EXPECT_EQ(par1.fluidModel, par2.fluidModel);
    EXPECT_EQ(par1.relaxModel, par2.relaxModel);
    EXPECT_ARR_EQ(par1.density, par2.density, BIN_PHASE);
    EXPECT_ARR_EQ(par1.relaxTime, par2.relaxTime, BIN_PHASE);
    EXPECT_FLOAT_EQ(par1.latticeSpeed, par2.latticeSpeed);
    EXPECT_FLOAT_EQ(par1.alpha0, par2.alpha0);
    EXPECT_FLOAT_EQ(par1.ift, par2.ift);
    EXPECT_FLOAT_EQ(par1.thickness, par2.thickness);
    for (int i = 0; i < MAX_GRAYSCALE; ++i) {
        EXPECT_EQ(par1.grayscale[i].type, par2.grayscale[i].type);
        EXPECT_FLOAT_EQ(par1.grayscale[i].value, par2.grayscale[i].value);
        EXPECT_EQ(par1.grayscale[i].direction, par2.grayscale[i].direction);
        EXPECT_EQ(par1.grayscale[i].phase, par2.grayscale[i].phase);
    }
}

TEST_P(FileParameterTest, load_file){
    ASSERT_NE(LoadFile(GetParam().parameterFile, &fileContent), ERROR);
    EXPECT_STREQ(fileContent, GetParam().content);
}

TEST_P(FileParameterTest, sync_string){
    if (ID == 0) {
        fileContent = (char*) malloc((strlen(GetParam().content) + 1) * sizeof(char));
        strcpy(fileContent, GetParam().content);
    }
    ASSERT_NE(SyncString(&fileContent, 0), ERROR);
    EXPECT_STREQ(fileContent, GetParam().content);
}

TEST_P(FileParameterTest, parse_parameter){
    PARAMETER par = DEFAULT_PARAMETER;
    ASSERT_NE(ParseParameter(GetParam().content, &par), ERROR);
    CompareParameter(par, GetParam().par);
}

TEST_P(FileParameterTest, load_parameter){
    PARAMETER par = DEFAULT_PARAMETER;
    ASSERT_NE(LoadParameter(GetParam().parameterFile, &par), ERROR);
    CompareParameter(par, GetParam().par);
}

typedef FileParameterTest FileParameterDeathTest;

void EXPECT_STRING_OBJECT(cJSON* root, const char key[],
                          int(*parser)(const cJSON* model, PARAMETER* par), std::string value="")
{
    PARAMETER par;
    cJSON_AddNumberToObject(root, key, 0);
    EXPECT_EQ(parser(root, &par), ERROR);
    cJSON_DeleteItemFromObjectCaseSensitive(root, key);
    cJSON_AddStringToObject(root, key, "");
    if (value != ""){
        EXPECT_EQ(parser(root, &par), ERROR);
        cJSON_DeleteItemFromObjectCaseSensitive(root, key);
        cJSON_AddStringToObject(root, key, value.c_str());
    }
}

void EXPECT_NUMBER_OBJECT(cJSON* root, const char key[],
                          int(*parser)(const cJSON* model, PARAMETER* par))
{
    PARAMETER par;
    cJSON_AddStringToObject(root, key, "");
    EXPECT_EQ(parser(root, &par), ERROR);
    cJSON_DeleteItemFromObjectCaseSensitive(root, key);
    cJSON_AddNumberToObject(root, key, 0);
}

void EXPECT_NUMBER_ARRAY_OBJECT(cJSON* root, const char key[],
                                int(*parser)(const cJSON* model, PARAMETER* par),
                                DECIMAL value=0, size_t min_length=0)
{
    PARAMETER par;
    cJSON_AddStringToObject(root, key, "");
    EXPECT_EQ(parser(root, &par), ERROR);
    cJSON_DeleteItemFromObjectCaseSensitive(root, key);
    cJSON* array = cJSON_CreateArray();
    cJSON_AddItemToObject(root, key, array);
    if (min_length > 0){
        EXPECT_EQ(parser(root, &par), ERROR);
        for (size_t i = 0; i < min_length; ++i)
            cJSON_AddItemToArray(array, cJSON_CreateString(""));
    } else {
        cJSON_AddItemToArray(array, cJSON_CreateString(""));
        min_length = 1; // ensure at least 1 item in the array for testing.
    }
    EXPECT_EQ(parser(root, &par), ERROR);
    while (cJSON_GetArraySize(array) > 0)
        cJSON_DeleteItemFromArray(array, 0);
    for (size_t i = 0; i < min_length; ++i)
        cJSON_AddItemToArray(array, cJSON_CreateNumber(0));
    if (value != 0){
        EXPECT_EQ(parser(root, &par), ERROR);
        cJSON_DeleteItemFromArray(array, 0);
        cJSON_AddItemToArray(array, cJSON_CreateNumber(value));
    }
}

TEST_F(FileParameterDeathTest, parse_grid_parameter){
    cJSON* grid = cJSON_CreateArray();
    EXPECT_STRING_OBJECT(grid, "gridFile", ParseGridParameter);
    EXPECT_NUMBER_ARRAY_OBJECT(grid, "gridSize", ParseGridParameter, 1);
    EXPECT_NUMBER_OBJECT(grid, "nDirection", ParseGridParameter);
    cJSON_Delete(grid);
}

TEST_F(FileParameterDeathTest, parse_control_parameter){
    cJSON* control = cJSON_CreateObject();
    for (auto key : {"startStep", "endStep", "resultIntvl", "restartIntvl", "reportIntvl",
                     "maxError"})
        EXPECT_NUMBER_OBJECT(control, key, ParseControlParameter);
    cJSON_Delete(control);
}

TEST_F(FileParameterDeathTest, parse_output_parameter){
    cJSON* output = cJSON_CreateObject();
    EXPECT_STRING_OBJECT(output, "outPath", ParseOutputParameter);
    EXPECT_NUMBER_ARRAY_OBJECT(output, "reportROI", ParseOutputParameter);
    cJSON_Delete(output);
}

TEST_F(FileParameterDeathTest, parse_model_parameter){
    cJSON* model = cJSON_CreateObject();
    EXPECT_STRING_OBJECT(model, "fluidModel", ParseModelParameter, "COLOR");
    EXPECT_STRING_OBJECT(model, "relaxModel", ParseModelParameter, "SRT");
    EXPECT_NUMBER_ARRAY_OBJECT(model, "density", ParseModelParameter, 0, 2);
    EXPECT_NUMBER_ARRAY_OBJECT(model, "relaxTime", ParseModelParameter, 0, 2);
    EXPECT_NUMBER_OBJECT(model, "latticeSpeed", ParseModelParameter);
    EXPECT_NUMBER_OBJECT(model, "alpha0", ParseModelParameter);
    EXPECT_NUMBER_OBJECT(model, "ift", ParseModelParameter);
    EXPECT_NUMBER_OBJECT(model, "thickness", ParseModelParameter);
    cJSON_Delete(model);
}

TEST_F(FileParameterDeathTest, parse_grayscale_parameter){
    PARAMETER par;
    cJSON* grayscale = cJSON_CreateArray();

    cJSON* item = cJSON_CreateObject();
    cJSON_AddStringToObject(item, "index", "0");
    cJSON_AddItemToArray(grayscale, item);
    item = cJSON_CreateObject();
    cJSON_AddNumberToObject(item, "index", -1);
    cJSON_AddItemToArray(grayscale, item);
    item = cJSON_CreateObject();
    cJSON_AddNumberToObject(item, "index", 0);
    cJSON_AddNumberToObject(item, "type", 0);
    cJSON_AddItemToArray(grayscale, item);
    item = cJSON_CreateObject();
    cJSON_AddNumberToObject(item, "index", 0);
    cJSON_AddStringToObject(item, "type", "NA");
    cJSON_AddItemToArray(grayscale, item);
    item = cJSON_CreateObject();
    cJSON_AddNumberToObject(item, "index", 0);
    cJSON_AddStringToObject(item, "type", "VD");
    cJSON_AddItemToArray(grayscale, item);
    item = cJSON_CreateObject();
    cJSON_AddNumberToObject(item, "index", 0);
    cJSON_AddStringToObject(item, "type", "GV");
    cJSON_AddItemToArray(grayscale, item);
    item = cJSON_CreateObject();
    cJSON_AddNumberToObject(item, "index", 0);
    cJSON_AddStringToObject(item, "type", "GV");
    cJSON_AddNumberToObject(item, "value", 1);
    cJSON_AddItemToArray(grayscale, item);
    item = cJSON_CreateObject();
    cJSON_AddNumberToObject(item, "index", 0);
    cJSON_AddStringToObject(item, "type", "GV");
    cJSON_AddNumberToObject(item, "value", 1);
    cJSON_AddNumberToObject(item, "direction", 1);
    cJSON_AddItemToArray(grayscale, item);
    item = cJSON_CreateObject();
    cJSON_AddNumberToObject(item, "index", 0);
    cJSON_AddStringToObject(item, "type", "SD");
    cJSON_AddItemToArray(grayscale, item);

    while (cJSON_GetArraySize(grayscale) > 0) {
        EXPECT_EQ(ParseGrayscaleParameter(grayscale, &par), ERROR);
        cJSON_DeleteItemFromArray(grayscale, 0);
    }

    cJSON_Delete(grayscale);
}

TEST_F(FileParameterDeathTest, parse_parameter){
    PARAMETER par;
    EXPECT_EQ(ParseParameter("{", &par), ERROR);
    EXPECT_EQ(ParseParameter("{}", &par), ERROR);
#define GRID \
    "    \"grid\" : {\n" \
    "        \"gridFile\":     \"Grid20x12y.raw\",\n" \
    "        \"gridSize\":     [20, 12],\n" \
    "        \"nDirection\":   9\n" \
    "    }\n"
    EXPECT_EQ(ParseParameter("{" GRID "}", &par), ERROR);
#define CONTROL \
    "    \"control\" : {\n" \
    "        \"startStep\":    0,\n" \
    "        \"endStep\":      1000,\n" \
    "        \"resultIntvl\":  100,\n" \
    "        \"restartIntvl\": 500,\n" \
    "        \"reportIntvl\":  10,\n" \
    "        \"maxError\":     0.000001\n" \
    "    }\n"
    EXPECT_EQ(ParseParameter("{" GRID "," CONTROL "}", &par), ERROR);
#define OUTPUT \
    "    \"output\" : {\n" \
    "        \"outPath\":      \"./LBFlow-Result/Grid20x12y\",\n" \
    "        \"reportROI\":    [0, 98]\n" \
    "    }\n"
    EXPECT_EQ(ParseParameter("{" GRID "," CONTROL "," OUTPUT "}", &par), ERROR);
#define MODEL \
    "    \"model\" : {\n" \
    "        \"fluidModel\":   \"SCSP\",\n" \
    "        \"relaxModel\":   \"SRT\",\n" \
    "        \"density\":      [1],\n" \
    "        \"relaxTime\":    [1],\n" \
    "        \"latticeSpeed\": 1\n" \
    "    }\n"
    EXPECT_EQ(ParseParameter("{" GRID "," CONTROL "," OUTPUT "," MODEL "}", &par), ERROR);
}


TEST_F(FileParameterDeathTest, load_parameter){
    PARAMETER par;
    EXPECT_EQ(LoadParameter("file dose not exist.", &par), ERROR);
    EXPECT_EQ(LoadParameter(SOURCE_DIR"examples/Channel3D-24x24y100z.raw", &par), ERROR);
}

#endif // FILE_PARAMETER_TEST_H
