/*! \file mpifunc.h
 * \brief MPI functions.
 *
 * - MPIInit() to initialize the MPI execution envoriment.
 * - extern ID: the local rank.
 * - extern MPI_SIZE: the communication size.
 * - remember to call MPI_Finalize() to terminates the MPI envoriment.
 *
 * \author Yu Zhou
 */
#ifndef MPIFUNC_H
#define MPIFUNC_H
#include <stdint.h>
#include "mpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief MPIInit initailize the MPI execution envoriment.
 * \param argc - number of arguments passed to the program.
 * \param argv - arguments passed to the program.
 * \return MPI error value.
 */
int MPIInit(int *argc, char ***argv);

#ifdef __cplusplus
}
#endif

#endif // MPIFUNC_H
