/*! \file testmpifunc.h
 * \brief Setup the mpi envoriment for all tests.
 * \author Yu Zhou
 */
#ifndef TESTMPIFUNC_H
#define TESTMPIFUNC_H

#include "gtest/gtest.h"
#include "mpifunc.h"

extern "C" {
    extern int ID;
    extern int MPI_SIZE;
}

class MPIEnviroment : public ::testing::Environment{
private:
    int  *m_argc;
    char ***m_argv;
public:
    MPIEnviroment(int* argc, char ***argv)
        :m_argc(argc) ,m_argv(argv) {}
protected:
    virtual void SetUp() {
        ASSERT_FALSE(MPIInit(m_argc, m_argv));
    }
    virtual void TearDown() {
        ASSERT_FALSE(MPI_Finalize());
    }
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

#endif // TESTMPIFUNC_H
