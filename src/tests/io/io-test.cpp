/*! \file io-test.cpp */

#include "file-raw-test.h"
#include "file-parameter-test.h"
#include "../tests/parallel/mpi-wrapper-test.h"
#include "mpi-wrapper.h"

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    // TODO: remove this line
    ::testing::Environment* const mpi_env =
            AddGlobalTestEnvironment(new MPIEnviroment(&argc, &argv));
    return RUN_ALL_TESTS();
}
