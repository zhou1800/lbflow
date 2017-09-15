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
    int LoadFile(const char* filePathName, char* content, int maxFileLength);
    int SyncString(char str[], int mId);
    int ParseParameterString(const char str[], PARAMETER* par);
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
    "        \"gridSize\":     [\"20\", \"12\"],\n"
    "        \"nDirection\":   \"9\",\n"
    "        \"latticeLen\":   \"1e-6\"\n"
    "    },\n"
    "    \"control\" : {\n"
    "        \"startTime\":        \"0\",\n"
    "        \"endTime\":          \"1000\",\n"
    "        \"writeInterval\":    \"100\",\n"
    "        \"reportInterval\":   \"10\",\n"
    "        \"restartInterval\":  \"500\",\n"
    "        \"maxError\":         \"0.000001\"\n"
    "    },\n"
    "    \"output\" : {\n"
    "        \"outPath\":      \"./LBFlow-Result/Grid20x12y\",\n"
    "        \"reportDir\":    \"1\",\n"
    "        \"reportSlice\":  [\"1\", \"98\"]\n"
    "    },\n"
    "    \"model\" : {\n"
    "        \"fluidModel\":       \"SCSP\",\n"
    "        \"equilModel\":       \"SRT\",\n"
    "        \"rexlationTime\":    \"1\",\n"
    "        \"latticeSpeed\":     \"1\",\n"
    "        \"density\":          \"1\"\n"
    "    },\n"
    "    \"index\" : {\n"
    "        \"0\" : {\n"
    "            \"type\" :  \"VD\",\n"
    "            \"phase\" : \"0\"\n"
    "        },\n"
    "        \"50\" : {\n"
    "            \"type\" :  \"VB\",\n"
    "            \"value\" : \"0.00001\",\n"
    "            \"dir\" :   \"1\",\n"
    "            \"phase\" : \"0\"\n"
    "        },\n"
    "        \"75\" : {\n"
    "            \"type\" :  \"GV\",\n"
    "            \"value\" : \"0.001\",\n"
    "            \"dir\" :   \"1\",\n"
    "        },\n"
    "        \"100\" : {\n"
    "            \"type\" :  \"VD\",\n"
    "            \"phase\" : \"1\"\n"
    "        },\n"
    "        \"150\" : {\n"
    "            \"type\" :  \"VB\",\n"
    "            \"value\" : \"0.00001\",\n"
    "            \"dir\" :   \"1\",\n"
    "            \"phase\" : \"1\"\n"
    "        },\n"
    "        \"200\" : {\n"
    "            \"type\" :  \"PB\",\n"
    "            \"value\" : \"1\",\n"
    "            \"dir\" :   \"2\",\n"
    "            \"phase\" : \"0\"\n"
    "        },\n"
    "        \"250\" : {\n"
    "            \"type\" :  \"SD\",\n"
    "            \"value\" : \"90\",\n"
    "        },\n"
    "        \"255\" : {\n"
    "            \"type\" :  \"SD\",\n"
    "            \"value\" : \"160\",\n"
    "        }\n"
    "    }\n"
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
        "../../../../../../../DRLab-R/Grid20x12y",
        1,
        {1, 98},
        SCSP,
        SRT,
        {1, 0},
        1,
        0,
        {1, 0},
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
    "        \"dimension\":    [\"24\", \"24\", \"100\"],\n"
    "        \"directions\":   \"19\",\n"
    "        \"latticeLen\":   \"1e-6\"\n"
    "    },\n"
    "    \"control\" : {\n"
    "        \"startTime\":        \"0\",\n"
    "        \"endTime\":          \"1000\",\n"
    "        \"writeInterval\":    \"100\",\n"
    "        \"reportInterval\":   \"10\",\n"
    "        \"restartInterval\":  \"500\",\n"
    "        \"maxError\":         \"0.000001\"\n"
    "    },\n"
    "    \"output\" : {\n"
    "        \"outPath\":      \"./LBFlow-Result/Channel3D-24x24y100z\",\n"
    "        \"reportDir\":    \"1\",\n"
    "        \"reportSlice\":  [\"1\", \"98\"]\n"
    "    },\n"
    "    \"model\" : {\n"
    "        \"fluidModel\":           \"COLOR\",\n"
    "        \"equilModel\":           \"MRT\",\n"
    "        \"rexlationTime\":        [\"0.57\", \"1.215416667\"],\n"
    "        \"latticeSpeed\":         \"1\",\n"
    "        \"density\":              \"1\",\n"
    "        \"alpha0\":               \"0.33333333333333333333\",\n"
    "        \"density\":              [\"1\", \"0.831683168\"],\n"
    "        \"interfacialTension\":   \"0.085920139\",\n"
    "        \"interfaceThickness\":   \"0.7\"\n"
    "    },\n"
    "    \"index\" : {\n"
    "        \"0\" : {\n"
    "            \"type\":  \"VD\",\n"
    "            \"phase\": \"0\"\n"
    "        },\n"
    "        \"50\" : {\n"
    "            \"type\":  \"VB\",\n"
    "            \"value\": \"0.00001\",\n"
    "            \"dir\":   \"1\",\n"
    "            \"phase\": \"0\"\n"
    "        },\n"
    "        \"75\" : {\n"
    "            \"type\":  \"GV\",\n"
    "            \"value\": \"0.001\",\n"
    "            \"dir\":   \"1\",\n"
    "        },\n"
    "        \"100\" : {\n"
    "            \"type\":  \"VD\",\n"
    "            \"phase\": \"1\"\n"
    "        },\n"
    "        \"150\" : {\n"
    "            \"type\":  \"VB\",\n"
    "            \"value\": \"0.00001\",\n"
    "            \"dir\":   \"1\",\n"
    "            \"phase\": \"1\"\n"
    "        },\n"
    "        \"200\" : {\n"
    "            \"type\":  \"PB\",\n"
    "            \"value\": \"1\",\n"
    "            \"dir\":   \"2\",\n"
    "            \"phase\": \"0\"\n"
    "        },\n"
    "        \"250\" : {\n"
    "            \"type\":  \"SD\",\n"
    "            \"value\": \"90\",\n"
    "        },\n"
    "        \"255\" : {\n"
    "            \"type\":  \"SD\",\n"
    "            \"value\": \"160\",\n"
    "        }\n"
    "    }\n"
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
        "../../../../../../../DRLab-R/Bera-122x102y102z",
        1,
        {1, 98},
        COLOR,
        MRT,
        {0.57, 1.215416667},
        1,
        0.33333333333333333333,
        {1, 0.831683168},
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

TEST_P(TestFilePar, load_file){
    char fileContent[MAX_FILE_LEN];
    ASSERT_NE(LoadFile(GetParam().parameterFile, fileContent, MAX_FILE_LEN),
              ERROR);
    EXPECT_STREQ(fileContent, GetParam().content);
}

TEST_P(TestFilePar, sync_file_content){
    char fileContent[MAX_FILE_LEN];
    const char * randomText = "Here are some contents in the file";
    if (ID == 0)
        strcpy(fileContent, randomText);
    ASSERT_NE(SyncString(fileContent, 0), ERROR);
    EXPECT_STREQ(fileContent, randomText);
}

TEST_P(TestFilePar, parse_parameter_string_content){
    PARAMETER par;
    ASSERT_NE(ParseParameterString(GetParam().content, &par), ERROR);
}

TEST_P(TestFilePar, load_parameter){
    PARAMETER par;
    ASSERT_NE(LoadParameter(GetParam().parameterFile, &par), ERROR);
}

TEST(DeathTestFilePar, load_parameter){
    PARAMETER par;
    char* veryLongFile = SOURCE_DIR"Examples/Channel3D-24x24y100z.raw";
    if (ID == 0) {
        EXPECT_DEATH(LoadParameter("non existing file", &par),
                     "Can not open parameter file");
        EXPECT_DEATH(LoadParameter(veryLongFile, &par),
                     "File size exceeds maximum allowed size");
    }
}

#endif // TESTFILEPAR_H
