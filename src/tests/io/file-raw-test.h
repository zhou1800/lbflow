/*! \file file-raw-test.h */

#ifndef FILE_RAW_TEST_H
#define FILE_RAW_TEST_H

#include "file-raw.h"
#include "../testcfg.h"
#include "gtest/gtest.h"

class FileRawTest: public ::testing::Test{
protected:
    virtual void SetUp(){
        raw = NULL;
    }
    virtual void TearDown(){
        free(raw);
    }
    RAW* raw;
};

TEST_F(FileRawTest, read_raw_file)
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

    ASSERT_NE(LoadRawFile(rawFile, &raw, size), ERROR);
    for (int i = 0; i < size; ++i)
        EXPECT_EQ(raw[i], content[i]);

}

typedef FileRawTest FileRawDeathTest;

TEST_F(FileRawDeathTest, read_raw_file){
    EXPECT_EQ(LoadRawFile(SOURCE_DIR"examples/Grid2D", &raw, 1), ERROR);
    EXPECT_EQ(LoadRawFile(SOURCE_DIR"examples/Grid20x12y.raw", &raw, 8192), ERROR);
}

#endif // FILE_RAW_TEST_H
