/*! \file testmpifunc.h
 * \brief Setup the mpi envoriment for all tests.
 * \author Yu Zhou (zhouyu@live.com)
 */
#ifndef TESTMPIFUNC_H
#define TESTMPIFUNC_H

#include "gtest/gtest.h"

extern "C" {
    #include "mpifunc.h"
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

#endif // TESTMPIFUNC_H
