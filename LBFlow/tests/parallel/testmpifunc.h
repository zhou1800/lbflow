/*! \file testmpifunc.h
 * \author Yu Zhou
 */
#ifndef TESTMPIFUNC_H
#define TESTMPIFUNC_H

#include "gtest/gtest.h"

extern "C" {
    #include "mpifunc.c"
}

class TestMPIFunc: public ::testing::Test{
protected:
    virtual void SetUp(){
    }
    virtual void TearDown(){
    }
};

#endif // TESTMPIFUNC_H
