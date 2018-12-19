/* \file slice-test.h */

#ifndef SLICE_TEST_H
#define SLICE_TEST_H

#include "slice.h"
#include "gtest/gtest.h"

class SliceTest: public ::testing::Test{
protected:
    virtual void SetUp(){
    }
    virtual void TearDown(){
    }
};

TEST(SliceTest, index_to_location) {
    int gridSize2D[] = {4, 7};
    int location2D[2] = {0};
    EXPECT_EQ(IndexToLocation(10, gridSize2D, location2D), SUCCESS);
    int gridSize3D[] = {10, 3, 2};
    int location3D[3] = {0};
    EXPECT_EQ(IndexToLocation(10, gridSize3D, location3D), SUCCESS);
}

TEST(SliceTest, location_to_index) {
    int index = 0;
    int gridSize2D[] = {4, 7};
    int location2D[] = {2, 3};
    EXPECT_EQ(LocationToIndex(location2D, gridSize2D, &index), SUCCESS);
    int gridSize3D[] = {10, 3, 2};
    int location3D[] = {3, 2, 0};
    EXPECT_EQ(LocationToIndex(location3D, gridSize3D, &index), SUCCESS);
}

#endif // SLICE_TEST_H
