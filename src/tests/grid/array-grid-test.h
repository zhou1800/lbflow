/* \file array-grid-test.h */

#ifndef ARRAY_GRID_TEST_H
#define ARRAY_GRID_TEST_H

#include "array-grid.h"
#include "../gtest-helper.h"
#include "gtest/gtest.h"

class ArrayGrid: public ::testing::Test{
protected:
    virtual void SetUp(){
        indices = NULL;
    }
    virtual void TearDown(){
        free(indices);
    }
    int* indices;
};

TEST_F(ArrayGrid, slice_to_indices) {
    int nIndices;
    SLICE emptySlice {0, {0}, {0}};
    ASSERT_EQ(SliceToIndices(&emptySlice, &indices, &nIndices), ERROR);

    SLICE slice2D {2, {3, 4}, {10, 1}};
    int expected2DIndices[] {2, 3, 4, 5, 12, 13, 14, 15, 22, 23, 24, 25};
    ASSERT_EQ(SliceToIndices(&slice2D, &indices, &nIndices), SUCCESS);
    ASSERT_EQ(nIndices, 3 * 4);
    EXPECT_ARR_EQ(indices, expected2DIndices, nIndices);

    SLICE slice3D {1, {4, 3, 2}, {100, 10, 1}};
    int expected3DIndices[] {1, 2, 11, 12, 21, 22, 101, 102, 111, 112, 121, 122,
                           201, 202, 211, 212, 221, 222, 301, 302, 311, 312, 321, 322};
    ASSERT_EQ(SliceToIndices(&slice3D, &indices, &nIndices), SUCCESS);
    ASSERT_EQ(nIndices, 4 * 3 * 2);
    EXPECT_ARR_EQ(indices, expected3DIndices, nIndices);
}

TEST_F(ArrayGrid, slice_to_box) {
    SLICE slice = {0, {0}, {0}};
    BOX box;
    EXPECT_EQ(SliceToBox(&slice, &box), SUCCESS);
}

TEST_F(ArrayGrid, box_to_slice) {
    BOX box = {{0}, {0}};
    SLICE slice;
    EXPECT_EQ(BoxToSlice(&box, &slice), SUCCESS);
}

TEST_F(ArrayGrid, index_to_location) {
    int gridSize2D[] = {4, 7};
    int location2D[2] = {0};
    EXPECT_EQ(IndexToPosition(10, gridSize2D, location2D), SUCCESS);
    int gridSize3D[] = {10, 3, 2};
    int location3D[3] = {0};
    EXPECT_EQ(IndexToPosition(10, gridSize3D, location3D), SUCCESS);
}

TEST_F(ArrayGrid, location_to_index) {
    int index = 0;
    int gridSize2D[] = {4, 7};
    int location2D[] = {2, 3};
    EXPECT_EQ(PositionToIndex(location2D, gridSize2D, &index), SUCCESS);
    int gridSize3D[] = {10, 3, 2};
    int location3D[] = {3, 2, 0};
    EXPECT_EQ(PositionToIndex(location3D, gridSize3D, &index), SUCCESS);
}

#endif // ARRAY_GRID_TEST_H
