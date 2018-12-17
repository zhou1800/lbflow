/*! \file mpi-wrapper.c */

#include "mpi-wrapper.h"

int ID;         /*!< Local rank. */
int MPI_SIZE;   /*!< Communication size. */

int MPIInit(int* argc, char*** argv)
{
    int mpiError = MPI_Init(argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &ID);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_SIZE);
    return mpiError;
}
