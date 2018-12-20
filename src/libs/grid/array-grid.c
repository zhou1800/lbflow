/* \file array-grid.c */

#include "array-grid.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int SliceToIndices(const SLICE* slice, int** indices, int* nIndices) {
    ASSERT((*nIndices = slice->sizes[0]) != 0, "[SliceToIndices] Slice is empty.\n");
    int dimensions = 1;
    for (int i = 1; i < MAX_DIMENSION && slice->sizes[i] != 0; ++i) {
        *nIndices *= slice->sizes[i];
        ++dimensions;
    }
    *indices = (int*) malloc(*nIndices * sizeof(int));
    int* temp = (int*) malloc(dimensions * sizeof(int));
    memcpy(temp, slice->sizes, dimensions * sizeof(int));
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

int SliceToBox(const SLICE* slice, BOX* box) {
    return SUCCESS;
}

int BoxToSlice(const BOX* box, SLICE* slice) {
    return SUCCESS;
}

int IndexToPosition(int index, const int size[], int position[]) {
    return SUCCESS;
}

int PositionToIndex(const int position[], const int size[], int* index) {
    return SUCCESS;
}
