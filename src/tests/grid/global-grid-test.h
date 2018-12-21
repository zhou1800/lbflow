/* \file global-grid-test.h */

#ifndef GLOBAL_GRID_TEST_H
#define GLOBAL_GRID_TEST_H

#include "global-grid.h"
#include "gtest/gtest.h"

extern "C" {
    extern int ID;
    extern int MPI_SIZE;
    int Primes(int value, int primes[]);
    int FindMaximum(const float array[], int n);
    int SeparateCommunication(int processor, int image_size[], int proc_size[]);
    int EffectivePoints(int effective, RAW globalGrid[], SLICE localSlice);
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
    const size_t size2D[] = {20, 12};    // Y, X order
    const size_t size3D[] = {3, 8, 10};    // Z, Y, X order
    RAW* localGrid = NULL;
    SLICE* localSlices = NULL;
    if (ID == 0) {
        EXPECT_EQ(ScatterGrid(globalGrid, size2D, localGrid, localSlices), SUCCESS);
        EXPECT_EQ(ScatterGrid(globalGrid, size3D, localGrid, localSlices), SUCCESS);
    }
    else {
        EXPECT_EQ(ScatterGrid(NULL, NULL, localGrid, localSlices), SUCCESS);
        EXPECT_EQ(ScatterGrid(NULL, NULL, localGrid, localSlices), SUCCESS);
    }
}

TEST(GlobalGridTest, PrimeTest)
{
    int primes[10];
    EXPECT_EQ(3, Primes(12, primes));
    EXPECT_EQ(2, primes[0]);
    EXPECT_EQ(2, primes[1]);
    EXPECT_EQ(3, primes[2]);
    EXPECT_EQ(0, Primes(1, primes));
    EXPECT_EQ(2, Primes(25, primes));
    EXPECT_EQ(5, primes[0]);
    EXPECT_EQ(5, primes[1]);
    EXPECT_EQ(1, Primes(17, primes));
    EXPECT_EQ(17, primes[0]);
}

TEST(GlobalGridTest, FindMaximumTest)
{
    float array1[] = {1, 2, 3};
    float array2[] = {10, 2, 1};
    EXPECT_EQ(2, FindMaximum(array1, 3));
    EXPECT_EQ(0, FindMaximum(array2, 3));
}

TEST(GlobalGridTest, SeparateCommunicationTest)
{
    int proc[3];
    int image1[3] = {80, 90, 1};
    int image2[3] = {80, 90, 70};
    SeparateCommunication(1, image1, proc);
    EXPECT_EQ(1, proc[0]);
    EXPECT_EQ(1, proc[1]);
    EXPECT_EQ(1, proc[2]);
    SeparateCommunication(12, image1, proc);
    EXPECT_EQ(4, proc[0]);
    EXPECT_EQ(3, proc[1]);
    EXPECT_EQ(1, proc[2]);
    SeparateCommunication(17, image1, proc);
    EXPECT_EQ(1, proc[0]);
    EXPECT_EQ(17, proc[1]);
    EXPECT_EQ(1, proc[2]);
    SeparateCommunication(25, image1, proc);
    EXPECT_EQ(5, proc[0]);
    EXPECT_EQ(5, proc[1]);
    EXPECT_EQ(1, proc[2]);
    SeparateCommunication(32, image2, proc);
    EXPECT_EQ(4, proc[0]);
    EXPECT_EQ(4, proc[1]);
    EXPECT_EQ(2, proc[2]);
}

TEST(GlobalGridTest, EffectivePointsTest) {
    RAW globalGrid[] = {     // 240 lattices.
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
    SLICE localSlice = {
        0,
        {3, 8, 10},
        {80, 10, 1}
    };
    SLICE localSlice2 = {
        0,
        {1, 5, 5},
        {240, 20, 1}
    };
    EXPECT_EQ(4, EffectivePoints(50, globalGrid, localSlice));
    EXPECT_EQ(72, EffectivePoints(255, globalGrid, localSlice));
    EXPECT_EQ(8, EffectivePoints(100, globalGrid, localSlice2));
    EXPECT_EQ(5, EffectivePoints(255, globalGrid, localSlice2));
}

#endif // GLOBAL_GRID_TEST_H
