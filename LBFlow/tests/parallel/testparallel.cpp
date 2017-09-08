/*! \file testparallel.cpp
 * \author Yu Zhou (zhouyu@live.com)
 */
#include "testmpifunc.h"

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
