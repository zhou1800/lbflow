/*! \file testfileraw.h
 * \author Yu Zhou
 */
#ifndef TESTFILERAW_H
#define TESTFILERAW_H

#include "gtest/gtest.h"
#include "../testcfg.h"
#include "fileraw.h"

class TestFileRaw: public ::testing::Test{
protected:
    virtual void SetUp(){
    }
    virtual void TearDown(){
    }
};

TEST(TestFileRaw, read_raw_file)
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
        255,255,255,255,255,255,255,255,255,255,
        255,255,255,255,255,255,255,255,255,255
    };
    const int size = 20 * 12;

    RAW* raw = new RAW[size];
    ASSERT_NE(LoadRawfile(rawFile, raw, size), ERROR);
    for (int i = 0; i < size; ++i)
        EXPECT_EQ(raw[i], content[i]);
    delete raw;
}

typedef TestFileRaw DeathTestFileRaw;

TEST(DeathTestFileRaw, load_parameter){
    const char* wrongFile = SOURCE_DIR"examples/Grid2D";
    const char* shortFile = SOURCE_DIR"examples/Grid20x12y.raw";
    const int size = 8192;
    RAW* raw = new RAW[size];
    EXPECT_DEATH({LoadRawfile(wrongFile, raw, size); exit(ERROR);},
                 "Can not open raw file");
    EXPECT_DEATH({LoadRawfile(shortFile, raw, size); exit(ERROR);},
                 "Failed to read raw file");
    delete raw;
}

#endif // TESTFILERAW_H
