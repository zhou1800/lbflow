/* \file global-grid.c */

#include "global-grid.h"
#include "config.h"
#include <mpi.h>
#include <math.h>

extern int ID;          /* Processor id */
extern int MPI_SIZE;    /* Communication size. */

/*! 
 * \brief prime
 * \param[in] value
 * \param[out] primes
 * \return number of primes
*/
int Primes(int value, int primes[]) {
    int i = 2, j = 0;
    while (value > 1) {
        if (value % i ==0) {
            primes[j++] = i; 
            value /= i;
        } else {
            i++;
        }
    }
    return j;
}

/*! 
 * \brief prime
 * \param[in] array
 * \param[in] n - length of array
 * \return index of max item
*/
int FindMaximum(const float array[], int n) {
    int i, index = 0;
    for (i = 1; i < n; i++) {
        if (array[i] > array[index]) {
            index = i;
        }
    }
    return index;
}

/*!
 * \brief seperate communications into dimensions, follow the image shape
 * \param[in] value - number of communiations
 * \param[in] image_size[3] - shape of image, {z, y, x}
 * \param[out] proc_size[3] - shape of communications matrix, {z, y, x}
 * \return SUCCESS
 */
int SeparateCommunication(int value, int image_size[], int proc_size[]) {
    int i, max;
    int primes[64];
    float dim_size[MAX_DIMENSION];
    for (i = 0; i<MAX_DIMENSION; i++) {
        dim_size[i] = image_size[i];
        proc_size[i] = 1;
    }
    i = Primes(value, primes);
    while (i > 0) {
        max = FindMaximum(dim_size, MAX_DIMENSION);
        proc_size[max] *= primes[i-1];
        dim_size[max] /= primes[--i];
    }
    return SUCCESS;
}

/*! 
 * \brief calculate effective points
 * \param[in] effective - value of effective point
 * \param[in] globalGrid
 * \param[in] localSlice - the slice need to calculate effective points
 * \return the number of effective points
*/
int EffectivePoints(int effective, RAW globalGrid[], SLICE localSlice) {
    int i, j, k, result = 0, index = localSlice.start, start[2];
    for (i = 0; i < localSlice.sizes[0]; i++) {
        start[0] = index;
        for (j = 0; j < localSlice.sizes[1]; j++) {
            start[1] = index;
            for (k = 0; k < localSlice.sizes[2]; k++) {
                result += (globalGrid[index] == effective ? 1 : 0);
                index += localSlice.stride[2];
            }
            index = start[1] + localSlice.stride[1];
        }
        index = start[0] + localSlice.stride[0];
    }
    return result;
}

int ScatterGrid(const RAW* globalGrid, const int globalSize[], RAW* localGrid, SLICE* localSlices) {
    return SUCCESS;
}
