/*! \file config.h
 * \brief Global configurations.
 * \author Yu Zhou (zhouyu@live.com)
 */

/*! \mainpage Lattice Boltzman Fluid Flow Simulator (LBFlow)
 *
 * The LBFlow implement the Lattice Boltzmann Method (LBM).
 *
 * Installation
 * ============
 * Use cmake to compile.
 *
 * Libraries, Executables and Tests
 * ================================
 * Libraries
 * ---------
 * |Libraries  |Description                                                  |
 * |-----------|-------------------------------------------------------------|
 * |libParallel|Parallel computiong libraries.                               |
 * |libIO      |File input and output.                                       |
 *
 * Executables
 * -----------
 * |Executables|Description                                                 |
 * |-----------|------------------------------------------------------------|
 * |LBFlow     |The Lattice Boltzmann Simulator.                            |
 *
 * Tests
 * -----
 * |Tests       |Description                                      |
 * |------------|-------------------------------------------------|
 * |io          |Test the io library.                             |
 * |parallel    |Test the parallel library.                       |
 *
 * Third Party libraries
 * ---------------------
 * |Third Party libraries |Description                             |
 * |----------------------|----------------------------------------|
 * |[googletest]          |For test driven development.            |
 *
 * [googletest]:     https://github.com/google/googletest
 *
 * Languages & Tools
 * =================
 *  - The project is mainly written in c and [c++] (standard c++11).
 *  - [MPI] [OpenMP] and [CUDA] are used for parallel computing.
 *  - Use [cmake] to build the project.
 *  - Use [doxygen] with [dot] and [latex] to generate this document.
 *  - [Git] is used for version control.
 *
 * [c++]:     http://www.stroustrup.com/
 * [MPI]:     http://www.mpi-forum.org/
 * [OpenMP]:  http://www.openmp.org/
 * [CUDA]:    https://developer.nvidia.com/cuda-zone
 * [cmake]:   https://cmake.org/
 * [doxygen]: http://www.stack.nl/~dimitri/doxygen/
 * [dot]:     http://graphviz.org/
 * [latex]:   http://www.latex-project.org/
 * [Git]:     http://www.git-scm.com/
 */
#ifndef CONFIG_H
#define CONFIG_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#define DECIMAL float       /*!< Precision of the decimal value. */

/*!
 * \brief The MODEL_CONSTS enum defines constants for the simulator.
 */
enum MODEL_CONSTS{
    MAX_DIMENSION = 3,      /*!< Maximum number of dimensions. */
    NUMBER_OF_SLICE = 2,    /*!< Number of report slices. */
    BIN_PHASE = 2,          /*!< Number of binary phase. */
};

enum PROGRAM_STATUS{
    NOERR = 0,          /*!< Return value on successful function call. */
    ERROR = -1,         /*!< Return value on function call failure. */
};

/*!
 * If expr is FALSE, print an error message and return ERROR.
 * \param[in] expr - expression.
 * \param[in] ... - error message.
 */
#define VALIDATE(expr, ...)\
    if (!(expr)) { fprintf(stderr, __VA_ARGS__); return ERROR; }

#ifdef __cplusplus
}
#endif

#endif // CONFIG_H
