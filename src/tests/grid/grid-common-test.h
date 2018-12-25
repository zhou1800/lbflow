/* \file grid-common-test.h */

#ifndef GRID_COMMON_TEST_H
#define GRID_COMMON_TEST_H

#include "grid-common.h"
#include "config.h"
#include "gtest/gtest.h"

extern "C" {
    extern int D;
    extern int Q;
}

class GridCommon: public ::testing::Test{
protected:
    virtual void SetUp(){
    }
    virtual void TearDown(){
    }
};

TEST_F(GridCommon, init) {
    EXPECT_EQ(Init(MAX_DIMENSION + 1), ERROR);
    EXPECT_EQ(Init(3), SUCCESS);
    EXPECT_EQ(D, 3);
    EXPECT_EQ(Q, 19);
    EXPECT_EQ(Init(2), SUCCESS);
    EXPECT_EQ(D, 2);
    EXPECT_EQ(Q, 9);
}

#endif // GRID_COMMON_TEST_H
