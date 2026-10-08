/* \file grid-test.cpp */

#include "grid-common-test.h"
#include "array-grid-test.h"
#include "global-grid-test.h"
#include "../mpi-environment.h"

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    AddGlobalTestEnvironment(new MPIEnviroment(&argc, &argv));
    return RUN_ALL_TESTS();
}
