/*! \file gridspace.h
 * \author Yu Zhou (zhouyu@live.com)
 */
#ifndef GRIDSPACE_H
#define GRIDSPACE_H

#include "defs.h"

namespace DRL {

namespace GRID {

/*!
 * \brief The GridSpace class provides the fundamental definition of grid's
 *  dimention(D) and direction(Q) space.
 *
 * All the grids in the same D&Q space should be initialized with the same
 *  GridSpace object. The motivation for designing this class is to share common
 *  values between all grids in the same D&Q space without declearing static
 *  global variables, and to keep the flexibility of defining multiple D&Q space
 *  in a single program.
 *
 * Currently, only the following "DmQn" combinations are avaliable:
 *  D1Q3, D2Q5, D2Q9, D3Q19, and D3Q27. See image "fig.1" and "fig.2" for
 *  examples of D2Q9 and D3Q19.
 * \image html 2Dmodel.png "fig.1 D2Q9 scheme"
 * \image latex 2Dmodel.png "fig.1 D2Q9 scheme"
 * \image html 3Dmodel.png "fig.2 D3Q19 scheme"
 * \image latex 3Dmodel.png "fig.2 D3Q19 scheme"
 */
class GridSpace
{
public:
    /// \name Grid dimension and direction definition.
    ///@{
    const size_t D;        //!< Number of dimensions.
    const size_t Q;        //!< Number of directions.

    ///@}

private:
    /*!
     * \brief The QAttrGroup struct groups all direction's attributes.
     */
    struct QAttrGroup{
        vava_int_t i;
        vava_val_t p;
        va_val_t     norm;
        va_idx_t     op;
        va_idx_t     half;
        vava_idx_t adj;
        vava_idx_t rel;
        vava_idx_t reld;
        vava_idx_t perp;
        va_idx_t     axis;
        va_int_t     sign;
    } qAttr;    /*!<
                 * \brief A helper variable for initializing direction's
                 *  attributes once for all, and use the values to initialize
                 *  public const variables.
                 */


public:
    /// \name Grid direction attributes
    ///@{
    const vava_int_t& Qi;   /*!<
                                * \brief Directions' projection on each
                                *  dimension.
                                *
                                * \note qi must be 0, -1, and 1 to act as
                                *  switch in multiply. Therefore,
                                *  schemes like D1Q5 are not supported.
                                */
    const vava_val_t& Qp;   /*!<
                                * \brief Length of directions' projection on
                                *  each dimension.
                                *
                                * Assume direction length is unit length.
                                */
    const va_val_t&     Qnorm;//!< Euclidean norm of each direction.
    const va_idx_t&     Qop;  //!< Opposit direction of each direction.
    const va_idx_t&     Qhalf;//!< Half of directions (exclude 0).
    const vava_idx_t& Qadj; //!< Adjecent directions of each direction.
    const vava_idx_t& Qrel; //!< Releated directions of each axis direction.
    const vava_idx_t& QrelD;//!< Releated directions of each axis.
    const vava_idx_t& Qprep;/*!<
                                * \brief Perpendicular directions of each axis
                                *  direction.
                                */
    const va_idx_t&     Qaxis;/*!<
                                * \brief Axis of directions orthogonal to axis
                                *  (avoid use direction 0).
                                */
    const va_int_t&     Qsign;/*!<
                                * \brief Sign of directions orthogonal to axis
                                *  (avoid use direction 0).
                                */
    ///@}
public:
    /*!
     * \brief Constructor.
     *
     * \sa initQ(size_t d, size_t q) and initQAttr(size_t d, size_t q) as
     *  constructor helper functions.
     * \param d - number of dimensions.
     * \param q - number of directions.
     */
    GridSpace(size_t d, size_t q);

private:
    /// \name Helper functions for constructor initializer list.
    ///@{
    /*!
     * \brief Perform initialize check for number of directions.
     *
     * If Q = 0 or Q > all avaliable directions, maxium directions will be used.
     *
     * \param d - number of dimensions.
     * \param q - number of directions.
     * \return checked number of directions.
     */
    static size_t initQ(size_t d, size_t q);

    /*!
     * \brief Calculate direction's attributes for initialization.
     * \param d - number of dimensions.
     * \param q - number of directions.
     * \return grouped directions's attributes.
     */
    static QAttrGroup initQAttr(size_t d, size_t q);
    ///@}
};

} // namespace GRID

} // namespace DRL

#endif // GRIDSPACE_H
