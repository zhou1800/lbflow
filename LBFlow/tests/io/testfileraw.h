/*! \file testfileraw.h
 * \author Yu Zhou
 */
#ifndef TESTFILERAW_H
#define TESTFILERAW_H

#include "gtest/gtest.h"

extern "C" {
    #include "fileraw.c"
}

class TestFileRaw: public ::testing::Test{
protected:
    virtual void SetUp(){
    }
    virtual void TearDown(){
    }
};

TEST_F(TestFileRaw, read_raw_file)
{
}

#endif // TESTFILERAW_H
