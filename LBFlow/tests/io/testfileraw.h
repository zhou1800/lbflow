/*! \file testfileraw.h
 * \author Yu Zhou (zhouyu@live.com)
 */
#ifndef TESTFILERAW_H
#define TESTFILERAW_H

#include "gtest/gtest.h"
#include "../testcfg.h"

extern "C" {
    #include "fileraw.h"
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
