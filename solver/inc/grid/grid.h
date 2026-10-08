/*! \file grid.h
 * \author Yu Zhou (zhouyu@live.com)
 */
#ifndef GRID_H
#define GRID_H

#include "gridspace.h"
#include <numeric>

namespace DRL {

namespace GRID {

/*!
 * \brief The Grid class represents an uniform rectilinear grid (image data).
 *
 * Advanced features:
 *  - Adjacent cells' indexes.
 *
 * \note The GridEx class add addintional attributes (as advanced features) to
 *  accerlate the calculation. Constructing these attributes are not only TIME
 *  CONSUMING, but also MEMORY CONSUMING. To be more efficient and flexible,
 *  these attibutes are calculated on the fly. Reference should be used to
 *  avoid checking the status of these attributes repetitively in a loop.
 */
class Grid
{
private:
    GridSpace* pS;          //!< This pointer to GridSpace is created on demand.

public:
    const GridSpace& S;       //!< Reference to the grid space.

    /// \name Grid attributes
    ///@{
    const va_idx_t Sz;        //!< Size of each dimension, in {z, y, x} order.
    const va_idx_t Stride;    //!< Strides of the grid (in gslice style).
    const size_t   N;         //!< Number of elements in the grid.

    ///@}

private:
    /// \name Advanced grid attributes (initialize on the fly).
    ///@{
    mutable vava_idx_t* adjI;    /*!<
                                     * \brief Adjacent cell's index on each
                                     *  direction.
                                     *
                                     * \note For cache optimization,
                                     *  index = adjI[_i_][_q_].
                                     *
                                     *  Second Order            |First Order
                                     * -------------------------|--------------
                                     *  current cell's index [i]|direction [q]
                                     */
    ///@}

public:
    /// \name Constructors and destructor.
    ///@{
    /*!
     * \brief Construct with existing GridSpace.
     * \param s - GridSpace to be used.
     * \param sz - dimension size of the grid.
     */
    Grid(const GridSpace& s, const va_idx_t& sz);

    /*!
     * \brief Construct and create new GridSpace.
     * \param sz - dimension size of the grid.
     * \param q - number of directions.
     */
    Grid(const va_idx_t& sz, size_t q);

    /*!
     * \brief Destructor.
     */
    virtual ~Grid();
    ///@}

    /// \name Get attributes.
    ///@{
    /*!
     * \brief Get adjecent cell's index on each direction.
     *
     * The adjacent indexes are calculated on the fly, and cached after first
     *  use.
     *
     * \return adjecent cell's indexes in the order of [_i_][_q_].
     */
    const vava_idx_t& getAdjI() const;
    ///@}

public:
    /// \name Helper functions
    ///@{
    /*!
     * \brief Calculate index projections on each dimension.
     * \sa CalcP2I() for inverse operation.
     * \param i - origin index.
     * \param sz - size of grid.
     * \return array of projection indexes on each dimension.
     */
    static va_idx_t CalcI2P(size_t i, const va_idx_t& sz);

    /*!
     * \brief Calculate index based on projection indexes.
     * \sa CalcI2P() for inverse operation.
     * \param p - index's projections on each dimension.
     * \param stride - strides of the grid owns the index.
     * \return index.
     */
    static size_t CalcP2I(va_idx_t p, const va_idx_t& stride);

private:
    /*!
     * \brief Calcualte the stride of grid.
     * \param sz - grid size.
     * \return grid strides.
     */
    static va_idx_t calcStride(const va_idx_t& sz);

protected:
    /*!
     * \brief Calculate adjecent index.
     *
     * The adjacent index is looped, i.e. xEnd + 1 = 0, 0 - 1 = xEnd
     *
     * \param i - current index in grid.
     * \param q - direction.
     * \param g - grid owns the index.
     * \return adjecent index.
     */
    static size_t calcAdjI(size_t i, size_t q, const Grid& g);
    ///@}
};

inline
const vava_idx_t &Grid::getAdjI() const
{
    if (!adjI) { // init adjI
        adjI = new vava_idx_t(N);
        for (size_t i = 0; i < N; ++i) { // for each cell
            (*adjI)[i].resize(S.Q);
            for (size_t q = 0; q < S.Q; ++q)    // for each direction
                (*adjI)[i][q] = calcAdjI(i, q, *this);
        }
    }
    return *adjI;
}

inline
va_idx_t Grid::CalcI2P(size_t i, const va_idx_t &sz)
{
    va_idx_t ip(sz.size());
    for (int d = static_cast<int>(sz.size()) - 1; d >= 0; --d) {
        ip[d] = i % sz[d];
        i = (i - ip[d]) / sz[d];
    }
    return ip;
}

inline
size_t Grid::CalcP2I(va_idx_t p, const va_idx_t &stride)
{
    p *= stride;
    return std::accumulate(begin(p), end(p), static_cast<size_t>(0));
}

inline
va_idx_t Grid::calcStride(const va_idx_t &sz)
{
    va_idx_t stride(1, sz.size());
    for (size_t i = sz.size() - 1; i > 0; --i)
        stride[i-1] = stride[i] * sz[i];
    return stride;
}

inline
size_t Grid::calcAdjI(size_t i, size_t q, const Grid &g)
{
    size_t adjI{i};
    for (size_t d = 0; d < g.S.D; ++d){
        size_t blockSz = g.Stride[d] * g.Sz[d];
        adjI += (i + g.S.Qi[q][d] * g.Stride[d] + blockSz)
                % blockSz - i % blockSz;
    }
    return adjI;
}

} // namespace GRID

} // namespace DRL

#endif // GRID_H
