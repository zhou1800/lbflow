/* \file global-grid-test.h */

#ifndef GLOBAL_GRID_TEST_H
#define GLOBAL_GRID_TEST_H

#include "global-grid.h"
#include "gtest/gtest.h"

class GlobalGridTest: public ::testing::Test{
protected:
    virtual void SetUp(){
    }
    virtual void TearDown(){
    }
};

TEST(GlobalGridTest, MpiCommunicationSizeTest)
{
    EXPECT_EQ(1, GetCommunicationSize());
}

TEST(GlobalGridTest, SeparateCommunicationTest)
{
    int comm_x=0, comm_y=0, comm_z=0;
    SeparateCommunication(12, 90, 100, 1, &comm_x, &comm_y, &comm_z);
    EXPECT_EQ(3, comm_x);
    EXPECT_EQ(4, comm_y);
    EXPECT_EQ(1, comm_z);
    SeparateCommunication(12, 90, 100, 2, &comm_x, &comm_y, &comm_z);
    EXPECT_EQ(2, comm_x);
    EXPECT_EQ(3, comm_y);
    EXPECT_EQ(2, comm_z);
    SeparateCommunication(4, 90, 100, 2, &comm_x, &comm_y, &comm_z);
    EXPECT_EQ(2, comm_x);
    EXPECT_EQ(2, comm_y);
    EXPECT_EQ(1, comm_z);
    SeparateCommunication(1, 90, 100, 2, &comm_x, &comm_y, &comm_z);
    EXPECT_EQ(1, comm_x);
    EXPECT_EQ(1, comm_y);
    EXPECT_EQ(1, comm_z);
    SeparateCommunication(7, 90, 100, 2, &comm_x, &comm_y, &comm_z);
    EXPECT_EQ(1, comm_x);
    EXPECT_EQ(7, comm_y);
    EXPECT_EQ(1, comm_z);
}

#endif // GLOBAL_GRID_TEST_H
