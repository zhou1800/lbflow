# Lattice Boltzman Fluid Flow Simulator (LBFlow)

## Introduction

LBFlow is a high performance multi-phase flow simulator. 
It implements the Lattice Boltzmann Method (LBM).

## Installation

### Requirements

- Any [c++] complier in able to compile c and [c++] (standard c++11) code.
- [cmake] to build the project.
- [MPI] [OpenMP] and [CUDA] for parallel computing.

[c++]:     http://www.stroustrup.com/
[MPI]:     http://www.mpi-forum.org/
[OpenMP]:  http://www.openmp.org/
[CUDA]:    https://developer.nvidia.com/cuda-zone
[cmake]:   https://cmake.org/

### Optional
- [doxygen] with [dot] and [latex] to generate document.
- [Git] to get the latest version.

[doxygen]: http://www.stack.nl/~dimitri/doxygen/
[dot]:     http://graphviz.org/
[latex]:   http://www.latex-project.org/
[Git]:     http://www.git-scm.com/

## Libraries, Executables and Tests

### Executables

* **LBFlow** - The Lattice Boltzmann Simulator.

### Libraries

* **libParallel** - Parallel computiong libraries.
* **libIO** - File input and output.

### Third Party libraries

* **[googletest]** - For testing and test driven development.
* **[cJSON]** - Parse JSON format parameter file.

[googletest]:     https://github.com/google/googletest
[cJSON]:          https://github.com/DaveGamble/cJSON

### Tests

* **io** - Test the io library.
* **parallel** - Test the parallel library.

## Contact Information

Yu (Joey) Zhou
[private contact omitted]
