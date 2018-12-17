/*! \file parallel-test.cpp */

#include "mpi-wrapper-test.h"

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    // TODO: remove the following line.
    ::testing::Environment* const mpi_env =
            AddGlobalTestEnvironment(new MPIEnviroment(&argc, &argv));
    return RUN_ALL_TESTS();
}
