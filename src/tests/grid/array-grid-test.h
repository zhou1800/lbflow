/* \file array-grid-test.h */

#ifndef ARRAY_GRID_TEST_H
#define ARRAY_GRID_TEST_H

#include "array-grid.h"
#include "gtest/gtest.h"

class ArrayGrid: public ::testing::Test{
protected:
    virtual void SetUp(){
    }
    virtual void TearDown(){
    }
};

TEST(ArrayGrid, slice_to_indices) {
    SLICE slice = {0, {0}, {0}};
    int* indices;
    int nIndices;
    EXPECT_EQ(SliceToIndices(&slice, indices, &nIndices), SUCCESS);
}

TEST(ArrayGrid, index_to_location) {
    int gridSize2D[] = {4, 7};
    int location2D[2] = {0};
    EXPECT_EQ(IndexToPosition(10, gridSize2D, location2D), SUCCESS);
    int gridSize3D[] = {10, 3, 2};
    int location3D[3] = {0};
    EXPECT_EQ(IndexToPosition(10, gridSize3D, location3D), SUCCESS);
}

TEST(ArrayGrid, location_to_index) {
    int index = 0;
    int gridSize2D[] = {4, 7};
    int location2D[] = {2, 3};
    EXPECT_EQ(PositionToIndex(location2D, gridSize2D, &index), SUCCESS);
    int gridSize3D[] = {10, 3, 2};
    int location3D[] = {3, 2, 0};
    EXPECT_EQ(PositionToIndex(location3D, gridSize3D, &index), SUCCESS);
}

#endif // ARRAY_GRID_TEST_H
