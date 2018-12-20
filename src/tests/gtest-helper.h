/* gtest-helper.h */

#ifndef GTEST_HELPER_H
#define GTEST_HELPER_H

#include "gtest/gtest.h"

void EXPECT_ARR_EQ(const int a[], const int b[], int n){
    for (int i = 0; i < n; ++i)
        EXPECT_EQ(a[i], b[i]);
}

void EXPECT_ARR_EQ(const DECIMAL a[], const DECIMAL b[], int n){
    for (int i = 0; i < n; ++i)
        EXPECT_FLOAT_EQ(a[i], b[i]);
}

#endif // GTEST_HELPER_H
