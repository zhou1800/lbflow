/*! \file mpi-wrapper-test.h
 * \brief Setup the mpi envoriment for all tests.
 */

#ifndef MPI_WRAPPER_TEST_H
#define MPI_WRAPPER_TEST_H

#include "mpi-wrapper.h"
#include "gtest/gtest.h"

extern "C" {
    extern int ID;
    extern int MPI_SIZE;
}

// TODO: Separate MPIWrapperTest
class MPIEnviroment : public ::testing::Environment{
public:
    MPIEnviroment(int* argc, char ***argv)
        :_argc(argc) ,_argv(argv) {}
protected:
    virtual void SetUp() {
        ASSERT_EQ(MPIInit(_argc, _argv), 0);
    }
    virtual void TearDown() {
        ASSERT_EQ(MPI_Finalize(), 0);
    }
private:
    int  *_argc;
    char ***_argv;
};

TEST(MPIEnviroment, initialize_MPI) {
    ASSERT_GE(ID, 0);
    ASSERT_GT(MPI_SIZE, 0);
}

TEST(MPIEnviroment, bcast_a_message) {
    int a;
    const int fixed = 1984;
    if (ID == 0)
        a = fixed;
    MPI_Bcast(&a, 1, MPI_INT, 0, MPI_COMM_WORLD);
    EXPECT_EQ(a, fixed);
}

#endif // MPI_WRAPPER_TEST_H
