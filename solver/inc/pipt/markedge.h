/*! \file markedge.h
 * \author Yu Zhou
 */

#ifndef MARKEDGE_H
#define MARKEDGE_H

#include "imgproc.h"
#include <sstream>

namespace DRL {

namespace PIPT {

/*!
 * \brief The MarkEdge class marks the edge of given grid (parallel).
 *
 * Firstly initialize the multi-core grid, then load the raw file to be marked.
 *  Finally save the marked file and show the "Surface area" / "Volum".
 */
class MarkEdge : public ImgProc
{
private:
    /// \name Parameters for calculating S/V
    ///@{
    uint64_t sumEdge;       ///< Sum of pixles on edge.
    uint64_t sumSolid;      ///< sum of pixles on solid.

    ///@}

public:
    /// \name Constructors and destructor.
    ///@{
    /*!
     * \brief Construct with existing GridMc.
     * \param g - GrdMc to be used.
     */
    MarkEdge(const GridMc& g);

    /*!
     * \brief Construct and create new GridMc
     * \param sz - dimension size of the grid.
     * \param q - number of directions.
     * \param ghoSz - ghost layer size (default = 1).
     */
    MarkEdge(const va_idx_t& sz, size_t q, size_t ghoSz = 1);
    ///@}

    /// \name Main operation.
    ///@{
    /*!
     * \brief Mark edge and calculate number of pixles on edge and solid.
     *
     * Every value not equal solid or boundary is considered OTHER. The value of
     *  SOLID cell near the interface is changed to boundary.
     *
     * \sa non-const version of Mark(const valarray<T>& buff, const T &sd, const T
     *  &bd);
     * \param buff - raw file to be marked - will be modified.
     * \param sd - value of solid.
     * \param bd - value of boundary.
     */
    template<typename T>
    valarray<T>* Mark(valarray<T> *buff, const T &sd, const T &bd);

    /*!
     * \brief Mark edge and calculate number of pixles on edge and solid
     *  (const version)
     *
     * \sa const version of Mark(valarray<T> *buff, const T &sd, const T
     *  &bd);
     * \param buff - raw file to be marked - will be modified.
     * \param sd - value of solid.
     * \param bd - value of boundary.
     */
    template<typename T>
    valarray<T> Mark(valarray<T> buff, const T &sd, const T &bd)
    { return *Mark(&buff, sd, bd); }

    /*!
     * \brief Save buff to file.
     * \param buff - buffer to be saved.
     * \param fname - destination file name.
     * \return saved file name.
     */
    template<typename T>
    string SaveFile(valarray<T>* buff, const string& fname);

    ///@}

    /// \name Get functions.
    ///@{
    /*!
     * \brief Get sum of the edge pixles.
     * \return sum of the edge pixles.
     */
    size_t getSumEdge() const { return sumEdge; }

    /*!
     * \brief Get sum of the solid pixles.
     * \return sum of the solid pixles.
     */
    size_t getSumSolid() const { return sumSolid; }

    /*!
     * \brief Get the ratio of Surface pixles over volum pixles.
     * \return S/V.
     */
    double getSV() const
    { return (double)sumEdge / (double)(G.N - sumSolid); }
    ///@}

};

template<typename T>
valarray<T>* MarkEdge::Mark(valarray<T> *buff, const T &sd, const T &bd)
{
    uint64_t nEdge = 0;
    uint64_t nSolid = 0;
    const vava_idx_t& adj {G.LG.FG.getAdjI()};
    const va_idx_t& i2fMap {G.LG.getI2fMap()};
    for (size_t cur = 0; cur < G.LG.N; ++cur) {
        size_t curf = i2fMap[cur];
        if ((*buff)[curf] != sd) continue;
        nSolid++;
        for (size_t q = 1; q < G.S.Q; ++q){
            if ((*buff)[adj[curf][q]] != sd && (*buff)[adj[curf][q]] != bd){
                (*buff)[curf] = bd;
                nEdge++;
                break;
            }
        }
    }
    nSolid -= nEdge;
    mpi().Reduce(&nEdge, &sumEdge, 1, MPI_UINT64_T, MPI_SUM);
    mpi().Reduce(&nSolid, &sumSolid, 1, MPI_UINT64_T, MPI_SUM);

    return buff;
}

template <typename T>
std::string MarkEdge::SaveFile(valarray<T>* buff, const string &fname)
{
    char name[255];
    sprintf(name, "%s-EdgeMarked.raw", fname.c_str());
    FileRaw<T> rawFile(name);
    valarray<T> tmpGrid{(*buff)[G.LG.getI2fGs()]};
    rawFile.Save(G, tmpGrid);
    mpi().Barrier();

    return name;
}

} // namespace PIPT

} // namespace DRL

#endif // MARKEDGE_H
