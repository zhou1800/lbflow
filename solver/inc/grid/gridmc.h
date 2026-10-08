/*! \file gridmc.h
 * \author Yu Zhou (zhouyu@live.com)
 */
#ifndef MCGRID_H
#define MCGRID_H

#include "gridgho.h"
#include "thread.h"
#include <vector>

namespace DRL {

namespace GRID {

using std::vector;

/*!
 * \brief The GridMc class is grid for multi-core system.
 *
 * The GridMc class use GridGho to avoid conflicts.
 *
 * Features for multi-core manipulate:
 * |Feature                                                      |Global|Local|
 * |:------------------------------------------------------------|:----:|:---:|
 * |Core grid.                                                   |  Y   |     |
 * |Adjacent core id on each direction.                          |  O*  |  Y  |
 * |Index projections of core id on each dimension.              |  O   |  Y  |
 * |Grid (wgl) of local core.                                    |  O   |  Y  |
 * |Start point of local grid on global grid.                    |  O   |  Y  |
 * |Location of local grid on global grid as gslice.             |  O   |  Y  |
 * |Location of local grid (wgl) on global grid (wgl) as gslice. |  O   |  Y  |
 *
 * Features marked as O (optional) are advanced features. \n
 * * Optional adjacent core index is CG.getAdjI().
 *
 * \note Like the Grid class, advanced features cost huge amount of memory and
 *  time to construct (especially when complied with gcc). Therefor, attributes
 *  for these features are calculated on the fly.
 */
class GridMc : public GridGho
{
public:
    /// \name Core attributes
    ///@{
    Grid            CG;        //!< Core grid.
    const va_idx_t  LCp;       //!< Core index projection on each dimension.
    const va_idx_t  AdjCi;     //!< Adjecent core index on each direction.

    ///@}

    /// \name Grid attributes
    ///@{
    const gslice    LGs;       //!< Local grid on global grid as gslice.
    const gslice    LGsWgl;    //!< Local grid (wgl) distribution as gslice.
    const va_idx_t  LStart;    //!< Start position of local grid on glb grid.
    const GridGho   LG;        //!< Local grid.

    ///@}

private:

    /// \name Attributes for all cores (initialize on the fly).
    ///@{
    mutable
    vector<va_idx_t> subCp;    //!< Core indexes projection on each dimension.
    mutable
    vector<gslice>   subGs;    //!< Local grid distributions as gslices.
    mutable
    vector<gslice>   subGsWgl; //!< Local grid(wgl) distribution as gslice.
    mutable
    vector<va_idx_t> subStart; //!< Start positions of sub grids.
    mutable
    vector<GridGho>  subG;     //!< Local grids.

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
    GridMc(const GridSpace& gsp, const va_idx_t& sz, size_t ghoSz = 1);

    /*!
     * \brief Construct and create new GridSpace.
     * \param sz - dimension size of the grid.
     * \param q - number of directions.
     * \param ghoSz - ghost layer size (default = 1).
     */
    GridMc(const va_idx_t& sz, size_t q, size_t ghoSz = 1);

    /*!
     * \brief Destructor.
     */
    virtual ~GridMc();
    ///@}

    /// \name Get functions
    ///@{
    /*!
     * \brief Get core index's projections on each dimension.
     * \return array of core index's projections, i.e. {cz, cy, cx}.
     */
    const vector<va_idx_t>& getSubCp() const;

    /*!
     * \brief Get sub grid location as gslice.
     * \return array of sub grid location as gslice.
     */
    const vector<gslice>& getSubGs() const;

    /*!
     * \brief Get sub grid (wgl) location as gslice.
     * \return array of sub grid (wgl) location as gslice.
     */
    const vector<gslice>& getSubGsWgl() const;

    /*!
     * \brief Get sub grid start position in global grid.
     * \return array of sub grid start position.
     */
    const vector<va_idx_t>& getSubStart() const;

    /*!
     * \brief Get sub grid.
     * \return array of sub grid.
     */
    const vector<GridGho> &getSubG() const;
    ///@}

private:
    /// \name Helper functions for constructor initializer list.
    ///@{
    /*!
     * \brief Even core distribution on each dimension of the grid.
     *
     * Achieved by the following steps:
     *  1. Get the primitive numbers of total cores.
     *  2. Distribute the primitives to each dimension. Respecting the rule that
     *     the remained largest primitive always distributed to the dimension
     *     with longest segmentation.
     *
     * \exception std::invalid_argument on more cores than the elements in grid.
     * \return array of evend core distribution on each dimension.
     */
    va_idx_t initCoreDistribution() const;

    /*!
     * \brief Calculate adjacent core indexes on each direction.
     * \return array of adjacent core indexes.
     */
    va_idx_t initAdjCi() const;

    /*!
     * \brief Generate local grid gslice.
     *
     * This function return size balanced (unequal) local grides when the global
     *  grid could not be divided with no remainder by core size. It is also
     *  ensured that core 0 always get the smallest grid.
     *
     * \param cP - core index projections on each dimension.
     * \return gslice of local grid.
     */
    gslice initLclGs(const va_idx_t& cP) const;
    ///@}
};

inline
const vector<va_idx_t> &GridMc::getSubCp() const
{
    if (subCp.size() == 0) {
        for (int i = 0; i < mpi().getCnt(); ++i)
            subCp.push_back(CalcI2P(i, CG.Sz));
    }
    return subCp;
}

inline
const vector<gslice> &GridMc::getSubGs() const
{
    if (subGs.size() == 0) {
        const vector<va_idx_t>& cp {getSubCp()};
        for (int i = 0; i < mpi().getCnt(); ++i)
            subGs.push_back(initLclGs(cp[i]));
    }
    return subGs;
}

inline
const vector<gslice> &GridMc::getSubGsWgl() const
{
    if (subGsWgl.size() == 0) {
        const vector<gslice>& lclGs {getSubGs()};
        for (int i = 0; i < mpi().getCnt(); ++i)
            subGsWgl.push_back(addGho2Gs(lclGs[i], *this));
    }
    return subGsWgl;
}

inline
const vector<va_idx_t> &GridMc::getSubStart() const
{
    if (subStart.size() == 0) {
        const vector<gslice>& lclGs {getSubGs()};
        for (int i = 0; i < mpi().getCnt(); ++i)
            subStart.push_back(CalcI2P(lclGs[i].start(), Sz));
    }
    return subStart;
}

inline
const vector<GridGho> &GridMc::getSubG() const
{
    if (subG.size() == 0) {
        const vector<gslice>& subGs {getSubGs()};
        for (int i = 0; i < mpi().getCnt(); ++i)
            subG.push_back(GridGho(S, subGs[i].size(), GhoSz));
    }
    return subG;
}

} // namespace GRID

} //namespace DRL

#endif // MCGRID_H
