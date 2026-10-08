/*! \file roi.h
 * \author Yu Zhou
 */
#ifndef ROI_H
#define ROI_H

#include "gridgho.h"

namespace DRL {

namespace GRID {

using edge_gs_t = valarray<map<size_t, gslice>>;  //!< Grid's edge indexed by
                                                  //!<  [_qM_][_qS_] as gslice.

/*!
 * \brief The ROI class provides pre-defined region of interesting(ROI), and
 *  operations on these ROIs.
 *
 * The ROI are repersented as (slice), gslice, indrect array, and (mask array).
 *
 * Pre-defined ROI includes:
 *  - Gslices of outter boundary on each direction. (like an outter shadow)
 *  - Gslices of inner boundary on each direction. (like an inner shadow)
 *  - Gslices of shifted original grid on each direction.
 *  - Indirect arrays of grid cells with the same value.
 *  - Gslices of edge on each direction.
 * For more details, see functions discription.
 *
 * Operations of ROI include:
 *  - Intersection of gslice-gslice, and gslice-indirect array.
 */
class ROI
{
public:
    ROI();

    ///\name ROI operations
    ///@{
    /*!
     * \brief Calculate the intersection of an indirect array and a gslice.
     * \param idx - indirect array.
     * \param gs - gslice.
     * \param gsN - number of elements in the grid which the gslice belongs to.
     * \param pos - position of intersection elements in indirect array.
     * \return indirect array of the interesction.
     */
    static va_idx_t Intersect(va_idx_t& idx,
                               const gslice& gs, std::size_t gsN,
                               va_idx_t* pos = 0);

    /*!
     * \brief Calculate the intersection of two gslices.
     * \param gs1 - gslice 1.
     * \param gsN1 - number of elements in the grid which the gs1 belongs to.
     * \param gs2 - gslice 2.
     * \param gsN2 - number of elements in the grid which the gs2 belongs to.
     * \param pos - position of intersection elements in gslice 1.
     * \return indirect array of the interesction.
     */
    static va_idx_t Intersect(const gslice& gs1, size_t gsN1,
                               const gslice& gs2, size_t gsN2,
                               va_idx_t* pos = 0);

    /*!
     * \brief Calculate the intersection of two indirect array.
     * \param idx1- indirect array 1.
     * \param idx2- indirect array 2.
     * \param pos - position of intersection elements in indirection array 1.
     * \return indirect array of the interesction.
     */
    static va_idx_t Intersect(va_idx_t& idx1, va_idx_t& idx2,
                               va_idx_t* pos = 0);

    /*!
     * \brief Merge two indirect array.
     * \param idx1 - indirect array 1.
     * \param idx2 - indirect array 2.
     * \return merged indirect array.
     */
    static va_idx_t Merge(va_idx_t& idx1, va_idx_t& idx2);

    /*!
     * \brief Convert gslice to indirect array.
     * \param gs - the gslice
     * \param gsN - number of elements in the grid which the gslice belongs to.
     * \return valarray of an indirect array.
     */
    static va_idx_t gs2idx(const gslice& gs, size_t gsN);
    ///@}

    ///\name Get pre-defined RIO
    ///@{
    /*!
     * \brief Select inner boundary of a grid with ghost layer.
     *
     * Inner boundary chart (D2Q9 for example):
     *
     *     ooooooooo
     *     o6222225o
     *     o3000001o
     *     o3000001o
     *     o7444448o
     *     ooooooooo
     *
     * | direction | selection |
     * |----------:|:---------:|
     * |         0 | 012345678 |
     * |         1 |    158    |
     * |         2 |    625    |
     * |         3 |    367    |
     * |         4 |    748    |
     * |         5 |     5     |
     * |         6 |     6     |
     * |         7 |     7     |
     * |         8 |     8     |
     *
     * \param g - ghost grid.
     * \return array of gslices represents selections on each direction.
     */
    static va_gslice getIBd(const GridGho &g);

    /*!
     * \brief Select outter boundary of a grid with ghost layer.
     *
     * Outter boundary chart (D2Q9 for example):
     *
     *     622222225
     *     300000001
     *     300000001
     *     300000001
     *     300000001
     *     744444448
     *
     * | direction | selection |
     * |----------:|:---------:|
     * |         0 |     0     |
     * |         1 |     1     |
     * |         2 |     2     |
     * |         3 |     3     |
     * |         4 |     4     |
     * |         5 |     5     |
     * |         6 |     6     |
     * |         7 |     7     |
     * |         8 |     8     |
     *
     * \param g - ghost grid.
     * \return array of gslices represents selections on each direction.
     */
    static va_gslice getOBd(const GridGho &g);

    /*!
     * \brief Select shifted origional grid on each direction.
     *
     * shifted grid chart (direction 5 on D2Q9 for example):
     *
     *     OOsssssss
     *     Oxsssssss    O - ghost layer
     *     Oxsssssss    x - origin matrix
     *     Oxsssssss    s - shifted matrix
     *     OxxxxxxxO
     *     OOOOOOOOO
     *
     * \param g - ghost grid.
     * \return array of gslices represents selections on each direction.
     */
    static va_gslice getSft(const GridGho &g);

    /*!
     * \brief Get ROI by value of the raw file (look up table for cell index).
     *
     * Build indexes of the grid by values in the raw data.
     *
     * \param data - the raw data.
     * \return map of index value and grid indexes holds this same value.
     */
    template <typename T>
    static map_idx_t getRawIdx(const valarray<T>& data);

    /*!
     * \brief Get edge of the grid on each direction for shift.
     *
     * Grid's edge chart (D2Q9 for example):
     *
     *     6222225 (00  10     11  20     22  30     33  40     44 ) - Main&Sub
     *     3000001 (0-8 302647 158 407813 625 102548 367 206315 748) - Selection
     *     3000001 (50   51 52 55 60   62 63 66) - Main&Sub
     *     3000001 (7034 81 62 5  8041 52 73 6 ) - Selection
     *     3000001 (70   73 74 77 80   84 81 88) - Main&Sub
     *     7444448 (5012 63 84 7  6023 74 51 8 ) - Selection
     *
     * 1.Main direction defines the position of a divider. For example direction
     *   1 with {x, y} = {1, 0} cut the grid into 2 pieces in x direction by
     *   putting a divider between 158 and 024, direction 5 with {x, y} = {1, 1}
     *   cut the grid into 4 pieces by putting two dividers at 625-013 and
     *   158-024.
     * 2.Sub direction is to obtain certain parts in the divided grid. Take the
     *   grid divided main direction 5 for example, sub direction 1, 2, 5 select
     *   the region on direction 1(18), 2(26), 5(5), and sub direction 0 select
     *   the region with 0(0347).
     *
     * \param g - Grid to be manipulated.
     * \param sftSz - number of shifted layer on each direction.
     * \return valarray of mapped gslice, index of valarray for main direction,
     *  key of map for subDirection, gslice for selected grid edge.
     */
    static edge_gs_t getEdge(const Grid& g, const va_idx_t& sftSz);

    ///@}

private:
    /// \name Helper functions
    ///@{
    /*!
     * \brief Calculate start point of the selection based on markers.
     *
     * Marker divides the dimension by three start points, i.e. {negative,
     *  origin, and positive}. Direction q indicates which marker should be
     *  taken on each dimension. All markers from each dimensions together
     *  determine the start point of the selected region. Then the start index
     *  is caculated.
     *
     * \param q - direction.
     * \param marker - three start points on each dimension.
     * \return the start index of a gslice.
     */
    static size_t calcStart(const GridGho &g, size_t q,
                            const vava_idx_t &marker);

    /*!
     * \brief Calculate size of selection.
     *
     * Size of the selection depends on its direction. For non-zero direction's
     *  projection on dimension, the size is equal to ghose layer size. For zero
     *  direction's projection on dimension, the size is equal to the grid size
     *  on that direction.
     *
     * \param q - direction.
     * \return array of selection size.
     */
    static va_idx_t calcSize(const GridGho &g, size_t q);

    /*!
     * \brief Calculate the divider for dividing the grid into pieces by
     *  main direction q.
     * \param q - main direction.
     * \return array of divider.
     */
    static va_idx_t calcDivider(const Grid &g, const va_idx_t &sftSz,
                                 size_t q);

    /*!
     * \brief Calculate gslice of sub region.
     *
     * Select pieces from the grid. The grid is divided by divider. The
     *  selection is based on main direction and sub direction.
     * \param div - divider.
     * \param qM - main direction.
     * \param qS - sub direction.
     * \return gslice of region.
     */
    static gslice calcSlice(const Grid &g, const va_idx_t& div,
                            size_t qM, size_t qS);
    ///@}
};

template <typename T>
map_idx_t ROI::getRawIdx(const valarray<T> &data)
{
    map<size_t, vector<size_t>> tempBdy;
    for (size_t i = 0; i < data.size(); ++i) {
        tempBdy[data[i]].push_back(i);
    }

    map_idx_t idx;
    for (auto& a : tempBdy) {
        idx[a.first] = valarray<size_t>(&(a.second)[0], a.second.size());
    }
    return idx;
}

} // namespace GRID

} // namespace DRL

#endif // ROI_H
