/* \file global-grid-test.h */

#ifndef GLOBAL_GRID_TEST_H
#define GLOBAL_GRID_TEST_H

#include "global-grid.h"
#include "gtest/gtest.h"

class GlobalGridTest: public ::testing::Test{
protected:
    virtual void SetUp(){
    }
    virtual void TearDown(){
    }
};

TEST(GlobalGridTest, empty_test)
{
    ASSERT_TRUE(true);
}

#endif // GLOBAL_GRID_TEST_H
