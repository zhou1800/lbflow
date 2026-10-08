/*! \file defs.h
 *
 * Global definitions.
 *
 * \author Yu Zhou (zhouyu@live.com)
 */

/*! \mainpage Digital Rock Lab
 *
 * The Digital Rock Lab (DRLab) implement "Prallel Image Processing Toolbox
 *  (PIPT)" and "Lattice Boltzmann Fluid Simulatior (LBSim)".
 *
 * Installation
 * ============
 * Compile and run.
 *
 * Source file structure
 * =====================
 * Codes are organized in folders.
 * |Folder  |Contents                                                  |
 * |--------|----------------------------------------------------------|
 * |3rdparty|3rd party libraries.                                      |
 * |apps    |Applications.                                             |
 * |doc     |Document related files and doxygen project files.         |
 * |inc     |Include files.                                            |
 * |samples |Several sample data sets for demonstration and testing.   |
 * |src     |Source files.                                             |
 * |test    |Test codes.                                               |
 *
 * Libraries, Executables and Tests
 * ================================
 * Libraries
 * ---------
 * |Libraries |Description                                       |
 * |----------|--------------------------------------------------|
 * |algorithm |Provide mathmatical algothrims.                   |
 * |core      |Global definations shared in the project.         |
 * |file      |File input and output.                            |
 * |grid      |The grid data structure.                          |
 * |helper    |Helper functions.                                 |
 * |lbm       |The Lattice Boltzman Simulator.                   |
 * |parallel  |Wrapper of MPI interface.                         |
 * |pipt      |The Parallel Image Processing Toolbox.            |
 *
 * The grid data structure implements features of
 *  -# Basic grid and direction definations.
 *  -# Grid decomposition for parallel computing.
 *  -# Ghost layer added for data transfer in parallel computing.
 *  -# Region of interest (ROI) for grid segmentation and selection.
 *
 * Headers (.h) and source code(.cpp) for libraries are found saparately in
 *  inc/ and src/ folders.
 *
 * Executables
 * -----------
 * |Executables|Description                                      |
 * |-----------|-------------------------------------------------|
 * |Connected  |Calculate connect regions (under construction).  |
 * |Histogram  |Calculate the histogram of an image.             |
 * |LbSim      |The Lattice Boltzmann Simulator.                 |
 * |MarkEdge   |Mark the edge between pore-space and matrix.     |
 *
 * 3rd Party libraries
 * -------------------
 * |3rd Party libraries |Description                             |
 * |--------------------|----------------------------------------|
 * |Eigen               |Library for linear algebra.             |
 * |QuickLZ             |A fastest compression library.          |
 *
 * [Eigen]:     http://eigen.tuxfamily.org
 * [QuickLZ]:   http://www.quicklz.com/
 *
 * Tests
 * -----
 * |Tests      |Description                                      |
 * |-----------|-------------------------------------------------|
 * |tgrid      |Test the grid library.                           |
 * |tfile      |Test the file library.                           |
 *
 * Languages
 * =========
 *  - The language for coding this project is mainly [c++] (standard c++11).
 *  - [MPI] is used for parallel computing.
 *  - Use [cmake] to build the project.
 *  - Use [doxygen] with [dot] and [latex] to generate this document.
 *  - [GIT] is used for version control.
 *
 * [c++]:     http://www.stroustrup.com/
 * [MPI]:     http://www.mpi-forum.org/
 * [cmake]:   https://cmake.org/
 * [doxygen]: http://www.stack.nl/~dimitri/doxygen/
 * [dot]:     http://graphviz.org/
 * [latex]:   http://www.latex-project.org/
 * [git]:     http://www.git-scm.com/
 */

/*! \page abbr Notation:
 *
 * \brief List of abreviation and its alias in the project.
 *
 * Suffix:
 *  - i[0,1,2...] Scalar of Index.
 *  - x,y,z       Scalar of specific dimension.
 *  - p           Vector of projection on dimension.
 *  - n           Scalar of total number.
 *  - sz          Vector of size.
 *
 * Varibles:
 * - D    Grid dimension
 *    - D     Number of dimensions.   (Dn)
 * - Q    Direction.
 *    - Q     Number of direction.    (Qn)
 *    - Qi    Direction index.
 *    - Qx    Elements in direction X
 * - (e)  Elements / cell. (admitted.)
 *    - i     Element index. (in 1D array)
 *    - x     Element index on dimension X.
 *    - n     Number of elements.
 *    - np    Number of elements on dimension i.
 *    - nx    Number of elements on dimension X.
 *    - sz    Element size(on each dimension).
 * - c    Core / processor.
 *    - ci    Core index. (in 1D array)
 *    - cx    Core index on dimension X.
 *    - nc    Numbe of cores. (total)
 *    - ncp   Number of cores on dimension p.
 *    - ncx   Number of cores on dimension X.
 *    - csz   Core size (on each dimension).
 *
 * Others:
 * - wgl  For "with ghost layer"
 */

#ifndef DEFS_H
#define DEFS_H

#include "cfg.h"
#include <cstdint>
#include <valarray>
#include <vector>
#include <map>
#include <string>
#include <set>

namespace DRL {

using std::size_t;
using std::valarray;
using std::vector;
using std::gslice;
using std::map;
using std::set;
using std::string;

using va_idx_t   = valarray<size_t>;        //!< Array of indexes.
using vava_idx_t = valarray<va_idx_t>;      //!< Array of array of indexes.
using va_int_t   = valarray<int>;           //!< Array of integers.
using vava_int_t = valarray<va_int_t>;      //!< Array of array of integers.
using va_val_t   = valarray<double>;        //!< Array of value.
using vava_val_t = valarray<va_val_t>;      //!< Array of array of value.
using ve_idx_t   = vector<size_t>;          //!< Vector of indexes.
using veve_idx_t = vector<ve_idx_t>;        //!< Vector of vector of indexes.
using ve_val_t   = vector<double>;          //!< Vector of values.
using veve_val_t = vector<ve_val_t>;        //!< Vector of vector of values.
using map_idx_t  = map<size_t, va_idx_t>;   //!< Map of indexes.
using map_val_t  = map<size_t, va_val_t>;   //!< Map of values.
using va_gslice  = valarray<gslice>;        //!< Array of gslice.
using va_u8_t    = valarray<uint8_t>;       //!< Array of 8-bit unsigned data.
using va_u64_t   = valarray<uint64_t>;      //!< Array of 64-bit unsigned data.

} // namespace DRL

#endif // DEFS_H
