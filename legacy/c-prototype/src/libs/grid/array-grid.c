/* \file array-grid.c */

#include "array-grid.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int SliceToIndices(const SLICE* slice, size_t** indices, size_t* nIndices)
{
    ASSERT((*nIndices = slice->sizes[0]) != 0, "[SliceToIndices] Slice is empty.\n");
    size_t dimensions = 1;
    for (size_t i = 1; i < MAX_DIMENSION && slice->sizes[i] != 0; ++i) {
        *nIndices *= slice->sizes[i];
        ++dimensions;
    }
    *indices = (size_t*) malloc(*nIndices * sizeof(size_t));
    size_t* temp = (size_t*) malloc(dimensions * sizeof(size_t));
    memcpy(temp, slice->sizes, dimensions * sizeof(size_t));
    size_t start = slice->start;

    for (size_t i = 0; i < *nIndices; ++i) {
        (*indices)[i] = start;
        --temp[dimensions -1];
        start += slice->stride[dimensions - 1];
        for (size_t j = dimensions - 1; j && !temp[j]; --j) {
            start -= slice->stride[j] * slice->sizes[j];
            temp[j] = slice->sizes[j];
            --temp[j - 1];
            start += slice->stride[j - 1];
        }
    }
    free(temp);
    return SUCCESS;
}

int SliceToBox(const SLICE* slice, BOX* box)
{
    return SUCCESS;
}

int BoxToSlice(const BOX* box, SLICE* slice)
{
    return SUCCESS;
}

int IndexToPosition(size_t index, const size_t size[], size_t position[])
{
    return SUCCESS;
}

int PositionToIndex(const size_t position[], const size_t size[], size_t* index)
{
    return SUCCESS;
}
