/*! \file testfilepar.h
 * \author Yu Zhou
 */

#ifndef TESTFILEPAR_H
#define TESTFILEPAR_H

#include "gtest/gtest.h"
#include "../testcfg.h"
#include "filepar.h"
#include "mpifunc.h"

extern "C" {
    extern int ID;
    extern PARAMETER DEFAULT_PARAMETER;
    int LoadFile(const char* filePathName, char* content, int maxFileLength);
    int SyncString(char str[], int mId);
    int ParseParameter(const char str[], PARAMETER* par);
}

struct PARAMETER_FILE{
    const char* parameterFile;
    const char* content;
    PARAMETER par;
};

class TestFilePar : public ::testing::TestWithParam<PARAMETER_FILE>{
protected:
    virtual void SetUp(){
    }
    virtual void TearDown(){
    }
};

PARAMETER_FILE Grid2D{
    SOURCE_DIR"Examples/Grid2D.json",     // parameter file
    "{\n"
    "    \"grid\" : {\n"
    "        \"gridFile\":     \"Grid20x12y.raw\",\n"
    "        \"gridSize\":     [20, 12],\n"
    "        \"nDirection\":   9,\n"
    "        \"latticeLen\":   1e-6\n"
    "    },\n"
    "    \"control\" : {\n"
    "        \"startTime\":    0,\n"
    "        \"endTime\":      1000,\n"
    "        \"writeIntvl\":   100,\n"
    "        \"rptIntvl\":     10,\n"
    "        \"rstIntvl\":     500,\n"
    "        \"maxError\":     0.000001\n"
    "    },\n"
    "    \"output\" : {\n"
    "        \"outPath\":      \"./LBFlow-Result/Grid20x12y\",\n"
    "        \"rptDir\":       0\n"
    "    },\n"
    "    \"model\" : {\n"
    "        \"fluidModel\":   \"SCSP\",\n"
    "        \"equilModel\":   \"SRT\",\n"
    "        \"rlxTime\":      [1],\n"
    "        \"latticeSpd\":   1,\n"
    "        \"density\":      [1]\n"
    "    },\n"
    "    \"index\" : [\n"
    "        {\n"
    "            \"idx\":    0,\n"
    "            \"type\" :  \"VD\",\n"
    "            \"phase\" : 0\n"
    "        },\n"
    "        {\n"
    "            \"idx\":    50,\n"
    "            \"type\" :  \"VB\",\n"
    "            \"value\" : 0.00001,\n"
    "            \"dir\" :   1,\n"
    "            \"phase\" : 0\n"
    "        },\n"
    "        {\n"
    "            \"idx\":    75,\n"
    "            \"type\" :  \"GV\",\n"
    "            \"value\" : 0.001,\n"
    "            \"dir\" :   1,\n"
    "            \"phase\" : 0\n"
    "        },\n"
    "        {\n"
    "            \"idx\":    100,\n"
    "            \"type\" :  \"VD\",\n"
    "            \"phase\" : 1\n"
    "        },\n"
    "        {\n"
    "            \"idx\":    150,\n"
    "            \"type\" :  \"VB\",\n"
    "            \"value\" : 0.00001,\n"
    "            \"dir\" :   1,\n"
    "            \"phase\" : 1\n"
    "        },\n"
    "        {\n"
    "            \"idx\":    200,\n"
    "            \"type\" :  \"PB\",\n"
    "            \"value\" : 1,\n"
    "            \"dir\" :   2,\n"
    "            \"phase\" : 0\n"
    "        },\n"
    "        {\n"
    "            \"idx\":    250,\n"
    "            \"type\" :  \"SD\",\n"
    "            \"value\" : 90\n"
    "        },\n"
    "        {\n"
    "            \"idx\":    255,\n"
    "            \"type\" :  \"SD\",\n"
    "            \"value\" : 160\n"
    "        }\n"
    "    ]\n"
    "}\n",
    // parameters
    {
        "Grid20x12y.raw",
        {20, 12},
        2,
        9,
        1e-6,
        0,
        1000,
        100,
        10,
        500,
        0.000001,
        "./LBFlow-Result/Grid20x12y",
        0,
        {0, 0},
        SCSP,
        SRT,
        {1, 0},
        1,
        {1, 0},
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
    SOURCE_DIR"Examples/Channel3D.json",  // parameter file
    "{\n"
    "    \"grid\" : {\n"
    "        \"gridFile\":     \"Channel3D-24x24y100z.raw\",\n"
    "        \"gridSize\":     [24, 24, 100],\n"
    "        \"nDirection\":   19,\n"
    "        \"latticeLen\":   1e-6\n"
    "    },\n"
    "    \"control\" : {\n"
    "        \"startTime\":    0,\n"
    "        \"endTime\":      1000,\n"
    "        \"writeIntvl\":   100,\n"
    "        \"rptIntvl\":     10,\n"
    "        \"rstIntvl\":     500,\n"
    "        \"maxError\":     0.000001\n"
    "    },\n"
    "    \"output\" : {\n"
    "        \"outPath\":      \"./LBFlow-Result/Channel3D-24x24y100z\",\n"
    "        \"rptDir\":       1,\n"
    "        \"rptSlice\":     [1, 98]\n"
    "    },\n"
    "    \"model\" : {\n"
    "        \"fluidModel\":   \"COLOR\",\n"
    "        \"equilModel\":   \"MRT\",\n"
    "        \"rlxTime\":      [0.57, 1.215416667],\n"
    "        \"latticeSpd\":   1,\n"
    "        \"density\":      [1, 0.831683168],\n"
    "        \"alpha0\":       0.33333333333333333333,\n"
    "        \"ift\":          0.085920139,\n"
    "        \"intfThk\":      0.7\n"
    "    },\n"

    "    \"index\" : [\n"
    "        {\n"
    "            \"idx\":   0,\n"
    "            \"type\":  \"VD\",\n"
    "            \"phase\": 0\n"
    "        },\n"
    "        {\n"
    "            \"idx\":   50,\n"
    "            \"type\":  \"VB\",\n"
    "            \"value\": 0.00001,\n"
    "            \"dir\":   1,\n"
    "            \"phase\": 0\n"
    "        },\n"
    "        {\n"
    "            \"idx\":   75,\n"
    "            \"type\":  \"GV\",\n"
    "            \"value\": 0.001,\n"
    "            \"dir\":   1,\n"
    "            \"phase\": 0\n"
    "        },\n"
    "        {\n"
    "            \"idx\":   100,\n"
    "            \"type\":  \"VD\",\n"
    "            \"phase\": 1\n"
    "        },\n"
    "        {\n"
    "            \"idx\":   150,\n"
    "            \"type\":  \"VB\",\n"
    "            \"value\": 0.00001,\n"
    "            \"dir\":   1,\n"
    "            \"phase\": 1\n"
    "        },\n"
    "        {\n"
    "            \"idx\":   200,\n"
    "            \"type\":  \"PB\",\n"
    "            \"value\": 1,\n"
    "            \"dir\":   2,\n"
    "            \"phase\": 0\n"
    "        },\n"
    "        {\n"
    "            \"idx\":   250,\n"
    "            \"type\":  \"SD\",\n"
    "            \"value\": 90\n"
    "        },\n"
    "        {\n"
    "            \"idx\":   255,\n"
    "            \"type\":  \"SD\",\n"
    "            \"value\": 160\n"
    "        }\n"
    "    ]\n"
    "}\n",
    {
        "Channel3D-24x24y100z.raw",
        {24, 24, 100},
        3,
        19,
        1e-6,
        0,
        1000,
        100,
        10,
        500,
        0.000001,
        "./LBFlow-Result/Channel3D-24x24y100z",
        1,
        {1, 98},
        COLOR,
        MRT,
        {0.57, 1.215416667},
        1,
        {1, 0.831683168},
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

typedef TestFilePar DeathTestFilePar;

INSTANTIATE_TEST_CASE_P(2Dand3D, TestFilePar,
                        ::testing::Values(Grid2D, Channel3D));

void EXPECT_ARR_EQ(const int a[], const int b[], int n){
    for (int i = 0; i < n; ++i)
        EXPECT_EQ(a[i], b[i]);
}

void EXPECT_ARR_EQ(const DECIMAL a[], const DECIMAL b[], int n){
    for (int i = 0; i < n; ++i)
        EXPECT_FLOAT_EQ(a[i], b[i]);
}

void CompareParameter(const PARAMETER& par1, const PARAMETER& par2)
{
    EXPECT_STREQ(par1.gridFile, par2.gridFile);
    EXPECT_ARR_EQ(par1.size, par2.size, MAX_DIMENSION);
    EXPECT_EQ(par1.nDimension, par2.nDimension);
    EXPECT_EQ(par1.nDirection, par2.nDirection);
    EXPECT_FLOAT_EQ(par1.latticeLen, par2.latticeLen);
    EXPECT_EQ(par1.startTime, par2.startTime);
    EXPECT_EQ(par1.endTime, par2.endTime);
    EXPECT_EQ(par1.writeIntvl, par2.writeIntvl);
    EXPECT_EQ(par1.rptIntvl, par2.rptIntvl);
    EXPECT_EQ(par1.rstIntvl, par2.rstIntvl);
    EXPECT_FLOAT_EQ(par1.maxError, par2.maxError);
    EXPECT_STREQ(par1.outPath, par2.outPath);
    EXPECT_EQ(par1.rptDir, par2.rptDir);
    EXPECT_ARR_EQ(par1.rptSlice, par2.rptSlice, NUMBER_OF_SLICE);
    EXPECT_EQ(par1.fluidModel, par2.fluidModel);
    EXPECT_EQ(par1.equilModel, par2.equilModel);
    EXPECT_ARR_EQ(par1.rlxTime, par2.rlxTime, BIN_PHASE);
    EXPECT_FLOAT_EQ(par1.latticeSpd, par2.latticeSpd);
    EXPECT_ARR_EQ(par1.density, par2.density, BIN_PHASE);
    EXPECT_FLOAT_EQ(par1.alpha0, par2.alpha0);
    EXPECT_FLOAT_EQ(par1.ift, par2.ift);
    EXPECT_FLOAT_EQ(par1.intfThk, par2.intfThk);
    for (int i = 0; i < MAX_GRAY; ++i) {
        EXPECT_EQ(par1.index[i].type, par2.index[i].type);
        EXPECT_FLOAT_EQ(par1.index[i].value, par2.index[i].value);
        EXPECT_EQ(par1.index[i].dir, par2.index[i].dir);
        EXPECT_EQ(par1.index[i].phase, par2.index[i].phase);
    }
}

TEST_P(TestFilePar, load_file){
    char fileContent[MAX_FILE_LEN];
    ASSERT_NE(LoadFile(GetParam().parameterFile, fileContent, MAX_FILE_LEN),
              ERROR);
    EXPECT_STREQ(fileContent, GetParam().content);
}

TEST_P(TestFilePar, sync_file_content){
    char fileContent[MAX_FILE_LEN];
    if (ID == 0)
        strcpy(fileContent, GetParam().content);
    ASSERT_NE(SyncString(fileContent, 0), ERROR);
    EXPECT_STREQ(fileContent, GetParam().content);
}

TEST_P(TestFilePar, parse_parameter_string_content){
    PARAMETER par = DEFAULT_PARAMETER;
    ASSERT_NE(ParseParameter(GetParam().content, &par), ERROR);
    CompareParameter(par, GetParam().par);
}

TEST_P(TestFilePar, load_parameter){
    PARAMETER par = DEFAULT_PARAMETER;
    ASSERT_NE(LoadParameter(GetParam().parameterFile, &par), ERROR);
    CompareParameter(par, GetParam().par);
}

TEST(DeathTestFilePar, load_parameter){
    PARAMETER par;
    const char* veryLongFile = SOURCE_DIR"Examples/Channel3D-24x24y100z.raw";
    if (ID == 0) {
        EXPECT_DEATH(LoadParameter("non existing file", &par),
                     "Can not open parameter file");
        EXPECT_DEATH(LoadParameter(veryLongFile, &par),
                     "File size exceeds maximum allowed size");
    }
}

#endif // TESTFILEPAR_H
