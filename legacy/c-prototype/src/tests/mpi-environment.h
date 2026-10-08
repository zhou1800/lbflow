/*! \file mpi-environment.h
 * \brief Set up MPI environment for testing.
 */

#ifndef MPI_ENVIRONMENT_H
#define MPI_ENVIRONMENT_H

#include "mpi-wrapper.h"
#include "gtest/gtest.h"

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

#endif // MPI_ENVIRONMENT_H
