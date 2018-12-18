/* \file grid-test.cpp */

#include "global-grid-test.h"
#include "../mpi-environment.h"

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    AddGlobalTestEnvironment(new MPIEnviroment(&argc, &argv));
    return RUN_ALL_TESTS();
}
