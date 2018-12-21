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
    int* EffectivePoints(int effective, int dimension,
        const RAW* globalGrid, const SLICE* localSlice);
    SLICE* DivideEffectivePoints(int effective, int dimension, int parts,
        const RAW* globalGrid, const SLICE* localSlice);
    SLICE* DivideGridByDimension(int effective, int parts[],
        const RAW* globalGrid, const SLICE* localSlice);
}

class GlobalGridTest: public ::testing::Test{
protected:
    virtual void SetUp(){
    }
    virtual void TearDown(){
    }
};

TEST(GlobalGridTest, scatter_grid) {
    const RAW globalGrid[] = { // 240 lattices.
        //0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
        150,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,200, //  0
        150,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,200, //  1
        255,255,255,255,255,100,100,100,100,100,255,255,255,100,100,100,100,100,100,200, //  2
        150,100,100,100,100,100,100,100,100,255,255,255,255,255,100,100,100,100,100,200, //  3
        150,100,100,100,100,100,100,100,255,255,255,255,255,255,255,100,100,100,100,200, //  4
        255,255,255,255,255,100,100,100,255,255,255,255,255,255,255,100,100,100,100,200, //  5
        150,100,100,100,100,100,100,100,255,255,255,255,255,255,255,100,100,100,100,200, //  6
        150,100,100,100,100,100,100,100,100,255,255,255,255,255,100,100,100,100,100,200, //  7
        255,255,255,255,255,100,100,100,100,100,255,255,255,100,100,100,100,100,100,200, //  8
        150,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,200, //  9
        150,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,200, // 10
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255  // 11
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
    const RAW globalGrid[] = { // 240 lattices.
        //0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
        150,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,200, //  0
        150,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,200, //  1
        255,255,255,255,255,100,100,100,100,100,255,255,255,100,100,100,100,100,100,200, //  2
        150,100,100,100,100,100,100,100,100,255,255,255,255,255,100,100,100,100,100,200, //  3
        150,100,100,100,100,100,100,100,255,255,255,255,255,255,255,100,100,100,100,200, //  4
        255,255,255,255,255,100,100,100,255,255,255,255,255,255,255,100,100,100,100,200, //  5
        150,100,100,100,100,100,100,100,255,255,255,255,255,255,255,100,100,100,100,200, //  6
        150,100,100,100,100,100,100,100,100,255,255,255,255,255,100,100,100,100,100,200, //  7
        255,255,255,255,255,100,100,100,100,100,255,255,255,100,100,100,100,100,100,200, //  8
        150,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,200, //  9
        150,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,200, // 10
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255  // 11
    };
    const SLICE localSlice = {
        0,
        {1, 12, 20},
        {240, 20, 1}
    };
    int* result;
    result = EffectivePoints(100, 1, globalGrid, &localSlice);
    EXPECT_EQ(18, result[0]);
    EXPECT_EQ(13, result[3]);
    EXPECT_EQ(0, result[11]);
    FREE(result);
    result = EffectivePoints(255, 2, globalGrid, &localSlice);
    EXPECT_EQ(4, result[0]);
    EXPECT_EQ(6, result[9]);
    EXPECT_EQ(8, result[11]);
    EXPECT_EQ(1, result[19]);
    FREE(result);
}

TEST(GlobalGridTest, DivideEffectivePointsTest) {
    const RAW globalGrid[] = {     // 25 lattices.
        255,255,255,255,255,
        150,100,100,100,100,
        150,100,100,100,100,
        150,100,100,100,100,
        255,255,255,255,255,
    };
    const SLICE localSlice = {
        0,
        {1, 5, 5},
        {25, 5, 1}
    };
    const SLICE localSlice2 = {
        1,
        {4, 4, 1},
        {1, 5, 25}
    };
    SLICE* result;
    result = DivideEffectivePoints(100, 1, 3, globalGrid, &localSlice);
    EXPECT_EQ(0, result[0].start);
    EXPECT_EQ(2, result[0].sizes[1]);
    EXPECT_EQ(10, result[1].start);
    EXPECT_EQ(15, result[2].start);
    EXPECT_EQ(2, result[2].sizes[1]);
    FREE(result);
    result = DivideEffectivePoints(100, 0, 3, globalGrid, &localSlice2);
    EXPECT_EQ(1, result[0].start);
    EXPECT_EQ(2, result[0].sizes[0]);
    EXPECT_EQ(3, result[1].start);
    EXPECT_EQ(4, result[2].start);
    EXPECT_EQ(1, result[2].sizes[0]);
    FREE(result);
    result = DivideEffectivePoints(100, 1, 1, globalGrid, &localSlice);
    EXPECT_EQ(0, result[0].start);
    EXPECT_EQ(5, result[0].sizes[1]);
    FREE(result);
}

TEST(GlobalGridTest, DivideGridByDimensionTest)
{
    const RAW globalGrid[] = { // 240 lattices.
        //0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
        150,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,200, //  0
        150,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,200, //  1
        255,255,255,255,255,100,100,100,100,100,255,255,255,100,100,100,100,100,100,200, //  2
        150,100,100,100,100,100,100,100,100,255,255,255,255,255,100,100,100,100,100,200, //  3
        150,100,100,100,100,100,100,100,255,255,255,255,255,255,255,100,100,100,100,200, //  4
        255,255,255,255,255,100,100,100,255,255,255,255,255,255,255,100,100,100,100,200, //  5
        150,100,100,100,100,100,100,100,255,255,255,255,255,255,255,100,100,100,100,200, //  6
        150,100,100,100,100,100,100,100,100,255,255,255,255,255,100,100,100,100,100,200, //  7
        255,255,255,255,255,100,100,100,100,100,255,255,255,100,100,100,100,100,100,200, //  8
        150,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,200, //  9
        150,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,200, // 10
        255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255  // 11
    };
    const SLICE localSlice = {
        0,
        {1, 6, 20},
        {240, 20, 1}
    };
    int parts[3] = {1, 1, 1};
    int parts2[3] = {1, 2, 3};
    SLICE *result;
    result = DivideGridByDimension(0, parts, globalGrid, &localSlice);
    EXPECT_EQ(1, result[0].sizes[0]);
    EXPECT_EQ(6, result[0].sizes[1]);
    EXPECT_EQ(20, result[0].sizes[2]);
    FREE(result);
    result = DivideGridByDimension(100, parts2, globalGrid, &localSlice);
    EXPECT_EQ(0, result[0].start);
    EXPECT_EQ(8, result[1].start);
    EXPECT_EQ(14, result[2].start);
    EXPECT_EQ(60, result[3].start);
    EXPECT_EQ(66, result[4].start);
    EXPECT_EQ(76, result[5].start);
    FREE(result);
}

#endif // GLOBAL_GRID_TEST_H
