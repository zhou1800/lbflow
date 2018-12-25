/* \file global-grid.c */

#include "global-grid.h"
#include "config.h"
#include <mpi.h>
#include <math.h>
#include <stdlib.h>

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
    int i, max, prime;
    int primes[64];
    float dim_size[MAX_DIMENSION];
    for (i = 0; i<MAX_DIMENSION; i++) {
        dim_size[i] = image_size[i];
        proc_size[i] = 1;
    }
    i = Primes(value, primes);
    while (i > 0) {
        max = FindMaximum(dim_size, MAX_DIMENSION);
        prime = (dim_size[max] / primes[i-1] >= 1) ? primes[i-1] : dim_size[max];
        proc_size[max] *= prime;
        dim_size[max] /= prime;
        i--;
    }
    return SUCCESS;
}

/*! 
 * \brief count effective points in dimension
 * \param[in] effective - value of effective point
 * \param[in] dimension - divide which dimesion
 * \param[in] globalGrid
 * \param[in] localSlice - the slice need to calculate effective points
 * \param[out] result - the array of histogram
 * \return SUCCESS
*/
int EffectivePoints(int effective, int dimension,
    const RAW* globalGrid, const SLICE* localSlice, size_t **ppresult) {
    size_t i, j, k;
    size_t max_i, max_j, max_k;
    size_t step_i, step_j, step_k;
    size_t start_i, start_j;
    size_t index = localSlice->start;
    size_t* result;
    // calloc result array
    result = *ppresult = (size_t *) calloc(localSlice->sizes[dimension], sizeof(size_t));
    // initialize loop
    max_i = localSlice->sizes[dimension];
    max_j = localSlice->sizes[(dimension + 1) % 3];
    max_k = localSlice->sizes[(dimension + 2) % 3];
    step_i = localSlice->stride[dimension];
    step_j = localSlice->stride[(dimension + 1) % 3];
    step_k = localSlice->stride[(dimension + 2) % 3];
    // calculate histogram of dimension
    for (i = 0; i < max_i; i++) {
        start_i = index;
        for (j = 0; j < max_j; j++) {
            start_j = index;
            for (k = 0; k < max_k; k++) {
                if (globalGrid[index] == effective) {
                    result[i]++;
                }
                index += step_k;
            }
            index = start_j + step_j;
        }
        index = start_i + step_i;
    }
    return SUCCESS;
}

/*!
 * \brief divide effective points into several slices by single dimension
 * \param[in] effective - value of effective point
 * \param[in] dimension - divide which dimesion
 * \param[in] parts - divide how many parts
 * \param[in] globalGrid
 * \param[in] localSlice - the slice need to calculate effective points
 * \param[out] result - part slices
 * \return SUCCESS
*/
int DivideEffectivePoints(int effective, int dimension, int parts,
    const RAW* globalGrid, const SLICE* localSlice, SLICE** ppresult) {
    size_t i, last;
    size_t step = 0, partial = 0, sum = 0;
    size_t *histogram;
    SLICE* result;
    // calloc result array
    result = *ppresult = (SLICE *) malloc(parts * sizeof(SLICE));
    // return whole slice directly
    if (parts == 1) {
        result[0] = *localSlice;
        return SUCCESS;
    }
    // calloc histogram array
    EffectivePoints(effective, dimension, globalGrid, localSlice, &histogram);
    // calculate sum of effective points
    for (i = 0; i < localSlice->sizes[dimension]; i++) {
        sum += histogram[i];
    }
    for (i = 0; i < parts; i++) {
        last = step;
        while (partial < ((i + 1) * sum / parts)) {
            partial += histogram[step++];
        }
        result[i] = *localSlice;
        if (i == (parts - 1)) {
            result[i].sizes[dimension] = localSlice->sizes[dimension] - last;
        } else {
            result[i].sizes[dimension] = step - last;
        }
        result[i].start = localSlice->start + last * localSlice->stride[dimension];
    }
    // free histogram array
    FREE(histogram);
    return SUCCESS;
}

/*!
 * \brief divide gird by three dimension
 * \param[in] effective - value of effective point
 * \param[in] parts - divide how many parts
 * \param[in] globalGrid
 * \param[in] localSlice - divide which slice
 * \param[out] result - slices of result
*/
int DivideGridByDimension(int effective, int parts[],
    const RAW* globalGrid, const SLICE* localSlice, SLICE **ppresult) {
    int i, j, k, n;
    SLICE *first, *second, *third, *result;
    result = *ppresult = (SLICE *) calloc (parts[0] * parts[1] * parts[2], sizeof(SLICE));
    DivideEffectivePoints(effective, 0, parts[0], globalGrid, localSlice, &first);
    for (i = 0; i < parts[0]; i++) {
        DivideEffectivePoints(effective, 1, parts[1], globalGrid, &first[i], &second);
        for (j = 0; j < parts[1]; j++) {
            DivideEffectivePoints(effective, 2, parts[2], globalGrid, &second[j], &third);
            for (k = 0; k < parts[2]; k++) {
                n = (i * parts[1] + j) * parts[2] + k;
                result[n] = third[k];
            }
            FREE(third);
        }
        FREE(second);
    }
    FREE(first);
    return SUCCESS;
}

int ScatterGrid(const RAW* globalGrid, const size_t globalSize[],
    RAW* localGrid, SLICE* localSlices) {
    return SUCCESS;
}
