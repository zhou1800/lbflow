/*! \file gridgho.h
 * \author Yu Zhou
 */
#ifndef GHOGRID_H
#define GHOGRID_H

#include "grid.h"

namespace DRL {

namespace GRID {

/*!
 * \brief The GridGho class is a Grid with ghost layer(wgl).
 *
 * Advanced Features for manipulating grid(wgl):
 *  - Access to inner grid (inherited) and full grid (as a variable)
 *  - Gslice for mapping from inner to full grid.
 *  - Map of inner grid indexes to full grid indexes.
 *
 * Ghost layers chart:
 *
 *     OOOOOOOOO
 *     OxxxxxxxO	O - ghost layer
 *     OxxxxxxxO 	x - inner grid
 *     OxxxxxxxO
 *     OxxxxxxxO
 *     OOOOOOOOO
 *
 * \note The full-grid(wgl) is stored as a variable for compatibility to
 *  non-ghost grid.
 *
 * \note Like the Grid class, Advanced features cost huge amount of memory and
 *  time to construct. Therefore, addintional values for these features are
 *  calculated on the fly.
 *
 * \bug __The gslice has inconsistence behavior with gcc and msvs.__
 *  While compiled with gcc (4.8.3-9), std::gslice requires as much memory as
 *  the selection array size, (which is huge). But it runs much faster than
 *  compiled with msvc(2015).
 */
class GridGho : public Grid
{
public:
    /// \name Ghost grid attributes.
    ///@{
    Grid         FG;        //!< Full grid(wgl).
    const size_t GhoSz;     //!< Size of ghost layer.

    ///@}

private:
    /// \name Mapping attributes (initialize on the fly).
    ///@{
    mutable gslice*    i2fGs;   //!< Inner grid gslice in full-grid.
    mutable va_idx_t* i2fMap;  //!< Map of inner grid index to full grid index.

    ///@}

public:
    /// \name Constructors and destructor.
    ///@{
    /*!
     * \brief Construct with existing GridSpace.
     * \param gsp - GridSpace to be used.
     * \param sz - dimension size of the grid.
     * \param ghoSz - ghost layer size (default = 1).
     */
    GridGho(const GridSpace& gsp, const va_idx_t& sz, size_t ghoSz = 1);

    /*!
     * \brief Construct and create new GridSpace.
     * \param sz - dimension size of the grid.
     * \param q - number of directions.
     * \param ghoSz - ghost layer size (default = 1).
     */
    GridGho(const va_idx_t& sz, size_t q, size_t ghoSz = 1);

    /*!
     * \brief Destructor.
     */
    virtual ~GridGho();
    ///@}

    /// \name Get attributes.
    ///@{
    /*!
     * \brief Get gslice of inner-grid in full-grid.
     * \return gslice of inner-grid in full-grid.
     */
    const gslice& getI2fGs() const;

    /*!
     * \brief Get map of inner grid index to full-grid index.
     * \return array of index transfer map.
     */
    const va_idx_t& getI2fMap() const;
    ///@}

protected:
    /// \name Helper functions
    ///@{
    /*!
     * \brief Calculate the index from inner grid to full grid.
     * \param iIdx - inner grid index.
     * \param iSz - inner grid size.
     * \param gSz - ghost layer size.
     * \param fStride - full grid stride.
     * \return index in full grid.
     */
    static size_t calcI2f(size_t iIdx, const va_idx_t &iSz, size_t gSz,
                          const va_idx_t &fStride);

    /*!
     * \brief Add ghost layer to a gslice.
     *
     * Convert a gslice of original inner grid to a gslice in full grid.
     *
     * \param iGs - inner gslice
     * \param g - ghost grid to work with.
     * \return gslice of grid(wgl) in full grid.
     */
    static gslice addGho2Gs(const gslice &iGs, const GridGho& g);
    ///@}
};

inline
const gslice &GridGho::getI2fGs() const
{
    if (!i2fGs) // init i2fGs
        i2fGs = new gslice(calcI2f(0, Sz, GhoSz, FG.Stride), Sz, FG.Stride);
    return *i2fGs;
}

inline
const va_idx_t &GridGho::getI2fMap() const
{
    if (!i2fMap) { // init i2fMap
        i2fMap = new va_idx_t(N);
        for (size_t i = 0; i < N; ++i)
            (*i2fMap)[i] = calcI2f(i, Sz, GhoSz, FG.Stride);
    }
    return *i2fMap;
}

inline
size_t GridGho::calcI2f(size_t iIdx, const va_idx_t& iSz, size_t gSz,
                        const va_idx_t& fStride)
{
    va_idx_t tmp {CalcI2P(iIdx, iSz)};
    tmp += gSz;
    tmp *= fStride;
    return std::accumulate(begin(tmp), end(tmp), static_cast<size_t>(0));
}

inline
gslice GridGho::addGho2Gs(const gslice &iGs, const GridGho &g)
{
    static size_t offset {std::accumulate(begin(g.FG.Stride), end(g.FG.Stride),
                                          static_cast<size_t>(0)) * g.GhoSz};
    size_t start {calcI2f(iGs.start(), g.Sz, g.GhoSz, g.FG.Stride) - offset};
    return gslice(start, iGs.size() + 2 * g.GhoSz, g.FG.Stride);
}

} // namespace GRID

} // namespace DRL

#endif // GHOGRID_H
