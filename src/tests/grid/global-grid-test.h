/* \file global-grid-test.h */

#ifndef GLOBAL_GRID_TEST_H
#define GLOBAL_GRID_TEST_H

#include "global-grid.h"
#include "gtest/gtest.h"

extern "C" {
    extern int ID;
    extern int MPI_SIZE;
    int SeparateCommunication(
        int comm_size,
        long image_size_x, long image_size_y, long image_size_z,
        int *comm_x, int *comm_y, int *comm_z);
}

class GlobalGridTest: public ::testing::Test{
protected:
    virtual void SetUp(){
    }
    virtual void TearDown(){
    }
};

TEST(GlobalGridTest, scatter_grid) {
    const RAW globalGrid[] = {     // 240 lattices.
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
    const int size2D[] = {20, 12};    // Y, X order
    const int size3D[] = {3, 8, 10};    // Z, Y, X order
    RAW* localGrid = NULL;
    SLICE* localSlices = NULL;
    int localSliceIndex = 0;
    EXPECT_EQ(ScatterGrid(globalGrid, size2D, localGrid, localSlices, &localSliceIndex), SUCCESS);
}

TEST(GlobalGridTest, SeparateCommunicationTest)
{
    int comm_x=0, comm_y=0, comm_z=0;
    SeparateCommunication(12, 90, 100, 1, &comm_x, &comm_y, &comm_z);
    EXPECT_EQ(3, comm_x);
    EXPECT_EQ(4, comm_y);
    EXPECT_EQ(1, comm_z);
    SeparateCommunication(12, 90, 100, 2, &comm_x, &comm_y, &comm_z);
    EXPECT_EQ(2, comm_x);
    EXPECT_EQ(3, comm_y);
    EXPECT_EQ(2, comm_z);
    SeparateCommunication(4, 90, 100, 2, &comm_x, &comm_y, &comm_z);
    EXPECT_EQ(2, comm_x);
    EXPECT_EQ(2, comm_y);
    EXPECT_EQ(1, comm_z);
    SeparateCommunication(1, 90, 100, 2, &comm_x, &comm_y, &comm_z);
    EXPECT_EQ(1, comm_x);
    EXPECT_EQ(1, comm_y);
    EXPECT_EQ(1, comm_z);
    SeparateCommunication(7, 90, 100, 2, &comm_x, &comm_y, &comm_z);
    EXPECT_EQ(1, comm_x);
    EXPECT_EQ(7, comm_y);
    EXPECT_EQ(1, comm_z);
}

#endif // GLOBAL_GRID_TEST_H
