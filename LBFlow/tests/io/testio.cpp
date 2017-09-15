/*! \file testio.cpp
 * \author Yu Zhou (zhouyu@live.com)
 */
#include "testfileraw.h"
#include "testfilepar.h"
#include "mpifunc.h"
#include "../tests/parallel/testmpifunc.h"

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::Environment* const mpi_env =
            AddGlobalTestEnvironment(new MPIEnviroment(&argc, &argv));
    return RUN_ALL_TESTS();
}
