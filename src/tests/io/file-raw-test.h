/*! \file file-raw-test.h */

#ifndef FILE_RAW_TEST_H
#define FILE_RAW_TEST_H

#include "file-raw.h"
#include "../testcfg.h"
#include "gtest/gtest.h"

class FileRawTest: public ::testing::Test{
protected:
    virtual void SetUp(){
    }
    virtual void TearDown(){
    }
};

TEST(FileRawTest, read_raw_file)
{
    const char* rawFile = SOURCE_DIR"examples/Grid20x12y.raw";
    const RAW content[] = {
        50,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,200,
        50,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,200,
        255,255,255,255,255,0,0,0,0,0,255,255,255,0,0,0,0,0,0,200,
        150,100,100,100,100,0,0,0,0,255,255,255,255,255,0,0,0,0,0,200,
        150,100,100,100,100,0,0,0,255,255,255,255,255,255,255,0,0,0,0,200,
        255,255,255,255,255,0,0,0,255,255,255,255,255,255,255,0,0,0,0,200,
        50,0,0,0,0,0,0,0,255,255,255,255,255,255,255,0,0,0,0,200,
        50,0,0,0,0,0,0,0,0,255,255,255,255,255,0,0,0,0,0,200,
        255,255,255,255,255,0,0,0,0,0,255,255,255,0,0,0,0,0,0,200,
        150,100,100,100,100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,200,
        150,100,100,100,100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,200,
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255
    };
    const int size = 20 * 12;

    RAW raw[size];
    ASSERT_NE(LoadRawfile(rawFile, raw, size), ERROR);
    for (int i = 0; i < size; ++i)
        EXPECT_EQ(raw[i], content[i]);
}

typedef FileRawTest FileRawDeathTest;

TEST(FileRawDeathTest, read_raw_file){
    const char* wrongFile = SOURCE_DIR"examples/Grid2D";
    const char* shortFile = SOURCE_DIR"examples/Grid20x12y.raw";
    const int size = 8192;
    RAW raw[size];
    EXPECT_DEATH({LoadRawfile(wrongFile, raw, size); exit(ERROR);},
                 "Can not open raw file");
    EXPECT_DEATH({LoadRawfile(shortFile, raw, size); exit(ERROR);},
                 "Failed to read raw file");
}

#endif // FILE_RAW_TEST_H
