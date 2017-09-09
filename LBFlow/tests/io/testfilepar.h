/*! \file testfilepar.h
 * \author Yu Zhou (zhouyu@live.com)
 */

#ifndef TESTFILEPAR_H
#define TESTFILEPAR_H

#include "gtest/gtest.h"
#include "../testcfg.h"
#include "filepar.h"

extern "C" {
    int LoadFile(const char* filePathName, char* content, int maxFileLength);
}

struct PARAMETER_FILE{
    const char* parameterFile;
    const char* content;
    int nLines;
    const char* lines[MAX_LINE_LEN];
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
    "",
    26,                                 // nlines
    // lines
    {
        "gridFile     \"Grid20x12y.raw\";",
        "dimension    20 12;",
        "directions   9;",
        "latticeLen   1e-6;",
        "startTime       0;",
        "endTime         1000;",
        "writeInterval   100;",
        "reportInterval  10;",
        "restartInterval 500;",
        "stopError       0.000001;",
        "outPath         \"../../../../../../../DRLab-R/Grid20x12y\";",
        "reportDir       1;",
        "reportSlice     1 98;",
        "fluidModel          \"scsp\";",
        "equilModel          \"SRT\";",
        "rexlationTime       1;",
        "latticeSpeed        1;",
        "density             1;",
        "0           \"VD\"    \"N/A\"       \"N/A\"   \"0\"",
        "50          \"VB\"    \"0.00001\"   \"1\"     \"0\"",
        "75          \"GV\"    \"0.001\"     \"1\"     \"N/A\"",
        "100         \"VD\"    \"N/A\"       \"N/A\"   \"1\"",
        "150         \"VB\"    \"0.00001\"   \"1\"     \"1\"",
        "200         \"PB\"    \"1\"         \"2\"     \"0\"",
        "250         \"SD\"    \"90\"        \"N/A\"   \"N/A\"",
        "255         \"SD\"    \"160\"       \"N/A\"   \"N/A\"",
    },
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
    "        \"gridFile\":     \"Channel3D-24x24y100z.raw\";\n"
    "        \"dimension\":    [\"24\", \"24\", \"100\"];\n"
    "        \"directions\":   \"19\";\n"
    "        \"latticeLen\":   \"1e-6\";\n"
    "    }\n"
    "}",
    34,                                 // nlines
    // lines
    {
        "gridFile     \"Channel3D-24x24y100z.raw\";",
        "dimension    24 24 100;",
        "directions   19;",
        "latticeLen   1e-6;",
        "startTime       0;",
        "endTime         1000;",
        "writeInterval   100;",
        "reportInterval  10;",
        "restartInterval 500;",
        "stopError       0.000001;",
        "outPath         \"../../../../../../../DRLab-R/Bera-122x102y102z\";",
        "reportDir       1;",
        "reportSlice     1 98;",
        "fluidModel          \"scsp\";",
        "equilModel          \"SRT\";",
        "rexlationTime       1;",
        "latticeSpeed        1;",
        "density             1;",
        "fluidModel          \"COLOR\";",
        "equilModel          \"MRT\";",
        "rexlationTime       0.57    1.215416667;",
        "latticeSpeed        1;",
        "alpha0              0.33333333333333333333;",
        "density             1       0.831683168;",
        "surfaceTension      0.085920139;",
        "interfaceThickness  0.7;",
        "0           \"VD\"    \"N/A\"       \"N/A\"   \"0\"",
        "50          \"VB\"    \"0.00001\"   \"1\"     \"0\"",
        "75          \"GV\"    \"0.001\"     \"1\"     \"N/A\"",
        "100         \"VD\"    \"N/A\"       \"N/A\"   \"1\"",
        "150         \"VB\"    \"0.00001\"   \"1\"     \"1\"",
        "200         \"PB\"    \"1\"         \"2\"     \"0\"",
        "250         \"SD\"    \"90\"        \"N/A\"   \"N/A\"",
        "255         \"SD\"    \"160\"       \"N/A\"   \"N/A\"",
    },
    // parameters
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
    ASSERT_NE(LoadFile(GetParam().parameterFile, fileContent, MAX_FILE_LEN), ERROR);
    EXPECT_STREQ(fileContent, GetParam().content);
}

TEST(DeathTestFilePar, load_file){
    char fileContent[MAX_FILE_LEN];
    EXPECT_DEATH({LoadFile("non existing file", fileContent, MAX_FILE_LEN); exit(-1);},
                 "Can not open parameter file");
    EXPECT_DEATH({LoadFile(Channel3D.parameterFile, fileContent, 1); exit(-1);},
                 "File size exceeds maximum allowed size");
}

TEST_P(TestFilePar, load_parameter){
    PARAMETER par;
    ASSERT_EQ(LoadParameter(GetParam().parameterFile, &par), NOERR);
}

#endif // TESTFILEPAR_H
