/*! \file mpi-wrapper-test.h
 * \brief Setup the mpi envoriment for all tests.
 */

#ifndef MPI_WRAPPER_TEST_H
#define MPI_WRAPPER_TEST_H

#include "mpi-wrapper.h"
#include "gtest/gtest.h"
#include <stdio.h>

extern "C" {
    extern int ID;
    extern int MPI_SIZE;
}

class MPIWrapperTest: public ::testing::Test{
protected:
    virtual void SetUp() {
    }
    virtual void TearDown() {
    }
};

TEST(MPIWrapperTest, initialize_MPI) {
    ASSERT_GE(ID, 0);
    ASSERT_GT(MPI_SIZE, 0);
}

TEST(MPIWrapperTest, bcast_a_message) {
    int a;
    const int fixed = 1984;
    if (ID == 0)
        a = fixed;
    MPI_Bcast(&a, 1, MPI_INT, 0, MPI_COMM_WORLD);
    printf("ID: %d, Size: %d, Value: %d\n", ID, MPI_SIZE, a);
    EXPECT_EQ(a, fixed);
}

#endif // MPI_WRAPPER_TEST_H
