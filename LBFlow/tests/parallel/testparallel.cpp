/*! \file testparallel.cpp
 * \author Yu Zhou
 */
#include "testmpifunc.h"

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
