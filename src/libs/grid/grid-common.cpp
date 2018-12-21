/* \file grid-common.cpp */

#include "grid-common.h"
#include "config.h"

int D = 0;      /*< Number of dimensions */
int Q = 0;      /*< Number of directions */

int Init(int dimension)
{
    ASSERT((unsigned int)dimension <= MAX_DIMENSION,
           "[Init] Number of dimensions %d exceed limit %d.\n", dimension, MAX_DIMENSION);
    D = dimension;
    switch (D) {
        case 2:
            Q = 9; break;
        case 3:
            Q = 19; break;
        default:
            Q = 0;
    }
    return SUCCESS;
}
