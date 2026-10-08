/*! \file imgproc.h
 * \author Yu Zhou (zhouyu@live.com)
 */

#ifndef IMGPROC_H
#define IMGPROC_H

#include "fileraw.h"
#include "gridmc.h"

namespace DRL {

namespace PIPT {

using namespace DRL::GRID;
using namespace DRL::IO;

class ImgProc
{
private:
    GridMc* pG;          //!< This pointer to GridMc is created on demand.

public:
    /// \name Grid defination.
    ///@{
    const GridMc& G;    ///< The global grid.

    ///@}
public:
    /// \name Constructors and destructor.
    ///@{
    /*!
     * \brief Construct with existing GridMc.
     * \param g - GrdMc to be used.
     */
    ImgProc(const GridMc& g);

    /*!
     * \brief Construct and create new GridMc
     * \param sz - dimension size of the grid.
     * \param q - number of directions.
     * \param ghoSz - ghost layer size (default = 1).
     */
    ImgProc(const va_idx_t& sz, size_t q, size_t ghoSz = 1);

    /*!
     * \brief Destructor.
     */
    virtual ~ImgProc();
    ///@}

    /// \name Image manipulation.
    ///@{
    /*!
     * \brief Load raw file.
     * \param fname - name of the file to be loaded.
     * \return loaded raw file.
     */
    template <typename T>
    valarray<T> LoadRaw(const string& fname);

    /*!
     * \brief Load raw file with ghost layer.
     * \param fname - name of the file to be loaded.
     * \param fill - value of extra boundary.
     * \return loaded raw file.
     */
    template <typename T>
    valarray<T> LoadRaw(const string &fname, const T &fill);
    ///@}
};

template <typename T> inline
valarray<T> ImgProc::LoadRaw(const string& fname)
{
    FileRaw<T> rawFile(fname);
    valarray<T> gridRaw = rawFile.Load(1, G);
    mpi().Barrier();
    return gridRaw;
}

template <typename T>
valarray<T> ImgProc::LoadRaw(const string& fname, const T& fill)
{
    FileRaw<T> rawFile(fname);
    valarray<T> gridRaw = rawFile.Load(G, fill, 0);
    mpi().Barrier();
    return gridRaw;
}

} // namespace PIPT

} // namespace DRL


#endif // IMGPROC_H
