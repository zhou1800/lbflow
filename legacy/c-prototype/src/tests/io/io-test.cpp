/*! \file io-test.cpp */

#include "file-raw-test.h"
#include "file-parameter-test.h"
#include "../mpi-environment.h"

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    AddGlobalTestEnvironment(new MPIEnviroment(&argc, &argv));
    return RUN_ALL_TESTS();
}
