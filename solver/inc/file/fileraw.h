/*! \file fileraw.h
 * \author Yu Zhou
 */
#ifndef FILERAW_H
#define FILERAW_H
#include "file.h"
#include "gridmc.h"

namespace DRL {

namespace IO {

using namespace DRL::GRID;
using std::runtime_error;
using std::slice;

/*!
 * \brief The FileRaw class is for loading and saving raw file.
 *
 *  - Use template to manipulate different file type, i.e. 8-bit, 16-bit.
 *  - Provide reloaded function Load() and cooresponding Save() to deal with
 *    both serialized and block data.
 */
template<typename T>
class FileRaw : public File
{
private:
    bool   isAppend;  //!< Flag of weather to use append mode in fwrite() and
                      //!< continue read from the end of last read.
    mutable size_t pos;       //!< Current position in fread().

public:
    /*!
     * \brief Constructor.
     * \param filename - file name of the raw file for input and output.
     * \param append - use append mode in fwrite() on true.
     */
    FileRaw(const string& filename, bool append = false);

    /// \name Single thread Load() and Save()
    ///@{
    /*!
     * \brief Load raw file (serialized).
     * \sa contrary operation: Save(const valarray<T> &buff) const.
     * \throw runtime_error on failed to open file.
     * \param n - number of elements in the file.
     * \return elements in the file as a valarray.
     */
    valarray<T> Load(size_t n) const;

    /*!
     * \brief Save buff to raw file (seralized).
     *
     * The single thread, serialized version of Save.
     *
     * \sa contrary operation: Load(size_t n) const.
     * \throw runtime_error on failed to open file.
     * \param buff - buffer to be saved
     */
    void Save(const valarray<T> &buff) const;
    ///@}

    /// \name Block Load() and Save()
    /*!
     * \brief Load raw file and scatter blocks to slave cores.
     *
     * Load file into a full grid in master core, then scatter to slave cores
     *  based on slave core gslice.
     *
     * \sa equivalent load: Load(size_t glbN, size_t lclN, const gslice& lclGs)
     *  const.
     * \sa contrary operation: Save(const MultiGrid& glbGrid, valarray<T>
     *  &lclBuff) const.
     * \param grid - descripter of global grid distribution on multi-cores.
     * \param tag - tag attached to the message.
     * \return scattered blocks for slave core.
     */
    valarray<T> Load(const GridMc &grid, int tag = MPI_ANY_TAG) const;

    /*!
     * \brief Load raw file and scatter blocks to slave cores by local gslice.
     *
     * Load file into a full grid in master procesor, then gather local gslice,
     *  distribute data to slave cores based on local gslice.
     *
     * \deprecated Prefer Load(const MultiGrid& glbGrid) const version which
     * scatter loaded blocks by global multi grid.
     *
     * \sa equivalent load: Load(const MultiGrid& glbGrid) const.
     * \sa contrary operation: Save(size_t glbN, valarray<T> &lclBuff, const
     *  gslice& lclGs) const.
     * \param glbN - global size of raw file.
     * \param lclN - grid size on slave core.
     * \param lclGs - local gslice on slave core.
     * \param tag - tag attached to the message.
     * \return scattered blocks for slave core.
     */
    valarray<T> Load(size_t glbN, size_t lclN, const gslice& lclGs,
                     int tag = MPI_ANY_TAG) const;

    /*!
     * \brief Load raw file and scatter blocks(wgl) to slave cores.
     *
     * Load the raw file into a full grid(wgl) (boundary is filled with value of
     *  ghoFill) on master core. Then scatter to slave core according to
     *  gslice(wgl) defined in global multi-ghost-grid.
     *
     * \param grid - global grid.
     * \param ghoFill - the value to be filled in ghost layer.
     * \param tag - tag attached to the message. Tag without default value for
     *  distinguishing with Load (GridMc &grid, int tag=MPI_ANY_TAG) const.
     * \return scattered blocks(wgl) for slave core.
     */
    valarray<T> Load(const GridMc& grid, const T& ghoFill, int tag) const;

    /*!
     * \brief Save local buffer (as block) into a global grid.
     *
     * Master core gather buffer from slave cores, then save to a file.
     *
     * \sa equivalent save: Save(size_t glbN, valarray<T> &lclBuff, const
     *  gslice& lclGs) const.
     * \sa contrary operation: Load(const MultiGrid& glbGrid) const.
     * \param g - global grid, provide grid structure data.
     * \param lclBuff - local buffer data to be saved.
     * \param tag - tag attached to the message.
     */
    void Save(const GridMc& g, valarray<T> &lclBuff,
              int tag = MPI_ANY_TAG) const;

    /*!
     * \brief Save local buffer (as block) into a global grid.
     *
     * Master core gather local gslice from slave, then gather local buffer.
     *  Save local buffer to global buffer according to local gslice. Finally
     *  save global buffer to file.
     *
     * \deprecated Prefer Save(const MultiGrid& glbGrid, valarray<T> &lclBuff)
     *  const for save local buffer as block.
     *
     * \sa equivalent save: Save(const MultiGrid& glbGrid, valarray<T> &lclBuff)
     *  const.
     * \sa contrary operation: Load(size_t glbN, size_t lclN, const gslice&
     *  lclGs) const.
     * \param glbN - size of global buffer.
     * \param[in] lclBuff - local buffer (data).
     * \param lclGs - local gslice.
     * \param tag - tag attached to the message.
     */
    void Save(size_t glbN, valarray<T> &lclBuff, const gslice& lclGs,
              int tag = MPI_ANY_TAG) const;
    ///@}

    /// \name Seralized Load() and Save()
    /*!
     * \brief Load raw file and serialise to slave cores (with length
     *  multiplier).
     *
     * Load file into master core. Then send serialized data to slave cores. The
     *  sent data length is the number defined in glbGrid multiple by the length
     *  multiplier.
     *
     * This function is specially designed to handle the f-function in Lattice
     *  Boltzmann model. F-function is several times larger than the grid, it
     *  needs to be serialized save and load for resrart purpose.
     *
     * \sa equivalent load: Load(size_t glbN, size_t lclN) const.
     * \param lenMp - length multiplier.
     * \param g - global grid of raw file.
     * \param tag - tag attached to the message.
     * \return seralized elements for slave core.
     */
    valarray<T> Load(size_t lenMp, const GridMc &g,
                     int tag = MPI_ANY_TAG) const;

    /*!
     * \brief Load raw file and serialise to slave cores.
     *
     * Load the data file into master core. Then master core wait for data size
     *  from slave core in sequence (from 1 to max). After received the size,
     *  send serialized data to slave core.
     *
     * \deprecated Prefer Load(const MultiGrid &glbGrid, size_t lenMul) for
     * scatter loaded data by global multi grid.
     *
     * \sa equivalent load: Load(size_t lenMp, const MultiGrid &glbGrid) const.
     * \sa contrary operation: Save(size_t glbN, valarray<T> &lclBuff) const.
     * \param glbN - size of full data.
     * \param lclN - size of serialized data on slave core.
     * \param tag - tag attached to the message.
     * \return serialized elements for slave core.
     */
    valarray<T> Load(size_t glbN, size_t lclN, int tag = MPI_ANY_TAG) const;

    /*!
     * \brief Save local buffer (seralized) into a file.
     *
     * Master core gather local buffer from slave in sequence (0 to max). Store
     *  local buffer into a global buffer, then save to file.
     *
     * \sa  contrary operation: Load(size_t glbN, size_t lclN) const
     * \param glbN - size of global grid.
     * \param lclBuff - local buffer to be saved.
     * \param tag - tag attached to the message.
     */
    void Save(size_t glbN, valarray<T> &lclBuff, int tag = MPI_ANY_TAG) const;

    /*!
     * \brief Save local buffer (seralized) with unknown total size into a file.
     *
     * Master core gather local buffer size from slave, then get the buffer.
     *  Store local buffer into a global buffer, then save to file.
     *
     * \sa  contrary operation: Load(size_t glbN, size_t lclN) const
     * \param lclBuff - local buffer to be saved.
     * \param tag - tag attached to the message. Tag with no dafault value for
     *  distinguishing with Save(const valarray<T> &buff) const.
     */
    void Save(valarray<T> &lclBuff, int tag) const;
    ///@}
};

template <typename T>
FileRaw<T>::FileRaw(const string &filename, bool append)
    :File(filename)
    , isAppend(append), pos(0)
{

}

template<typename T> inline
valarray<T> FileRaw<T>::Load(size_t n) const
{
    valarray<T> buff(n);

    FILE* file = fopen(fname.c_str(), "rb");
    if (file == 0)
        throw runtime_error("FileRaw::Load() - Can not open file for input.");

    if (isAppend)
        if (fseek(file, pos, SEEK_SET))
            throw runtime_error("FileRaw::Load() - End of file while seeking.");

    if (fread(&buff[0], sizeof(T), n, file) != n)
        throw runtime_error("FileRaw::Load() - Failed to read from raw file");

    if (isAppend) pos = ftell(file);

    fclose(file);
    return buff;
}

template<typename T>
void FileRaw<T>::Save(const valarray<T>& buff) const
{
    FILE* file;
    if (isAppend) file = fopen(fname.c_str(), "ab");
    else file = fopen(fname.c_str(), "wb");

    if (file == 0)
        throw runtime_error("FileRaw::Save() - Can not open file for output.");

    if (fwrite(&buff[0], sizeof(T), buff.size(), file) != buff.size())
        throw runtime_error("FileRaw::Save() - Failed to save to file");

    fclose(file);
}

template<typename T> inline
valarray<T> FileRaw<T>::Load(const GridMc &grid, int tag) const
{
    valarray<T> lclBuff(grid.LG.N);

    if (mpi().getId() == 0){
        // Load all contents into buffer in master core
        valarray<T> glbBuff = Load(grid.N);
        // Scatte to slave cores.
        const vector<gslice>& subGs {grid.getSubGs()};
        for (auto id = 1; id < mpi().getCnt(); ++id) {
            valarray<T> srcBuff = glbBuff[subGs[id]];
            mpi().Send(&srcBuff[0], srcBuff.size() * sizeof(T), MPI_BYTE, tag,
                    id);
        }
        lclBuff = glbBuff[grid.LGs];
    }
    else {
        // receive from master core
        mpi().Recv(&lclBuff[0], lclBuff.size() * sizeof(T), MPI_BYTE, tag, 0);
    }
    return lclBuff;
}

template<typename T> inline
valarray<T> FileRaw<T>::Load(size_t glbN, size_t lclN, const gslice &lclGs,
                             int tag) const
{
    valarray<T> lclBuff(lclN);     //loacal buffer

    GsliceOpen srcGso(lclGs);       // tempory gslice for mpi syncorize
    MPI_Datatype gsPacked;
    // prepare gslice package to be synced
    void* pFirst = mpi().Pack(&gsPacked, &srcGso);

    if (mpi().getId() == 0){
        // Load all contents into global buffer in master core
        valarray<T> glbBuff = Load(glbN);
        for (auto i = 1; i < mpi().getCnt(); ++i) {
            // Receive local gslice from slave cores.
            mpi().Recv(pFirst, 1, gsPacked, tag);
            int src = mpi().getSrc();
            gslice srcGs = srcGso.getGs();
            valarray<T> srcBuff = glbBuff[srcGs];
            mpi().Send(&srcBuff[0], srcBuff.size() * sizeof(T), MPI_BYTE, tag,
                    src);
        }
        lclBuff = glbBuff[lclGs];
    }
    else {
        mpi().Send(pFirst, 1, gsPacked, tag, 0);
        mpi().Recv(&lclBuff[0], lclN * sizeof(T), MPI_BYTE, tag, 0);
    }
    return lclBuff;
}

template<typename T> inline
valarray<T> FileRaw<T>::Load(const GridMc &grid, const T &ghoFill,
                             int tag) const
{
    valarray<T> lclBuff(1, grid.LG.FG.N);

    if (mpi().getId() == 0){
        // initialize global buffer filled with ghoFill
        valarray<T> glbBuff(ghoFill, grid.FG.N);
        // Load all contents into center of global buffer on master core
        glbBuff[grid.getI2fGs()] = Load(grid.N);
        const vector<gslice>& subGs {grid.getSubGsWgl()};
        for (auto id = 1; id < mpi().getCnt(); ++id) {
            valarray<T> srcBuff = glbBuff[subGs[id]];
            mpi().Send(&srcBuff[0], srcBuff.size() * sizeof(T), MPI_BYTE, tag,
                    id);
        }
        lclBuff = glbBuff[grid.LGsWgl];
    }
    else {
        mpi().Recv(&lclBuff[0], lclBuff.size() * sizeof(T), MPI_BYTE, tag, 0);
    }
    return lclBuff;
}

template<typename T> inline
void FileRaw<T>::Save(const GridMc &g, valarray<T>& lclBuff, int tag) const
{
    if (mpi().getId() == 0){
        valarray<T> glbBuff(g.N);      // global buffer
        glbBuff[g.LGs] = lclBuff;      // copy master local buffer

        const vector<gslice>& subGs {g.getSubGs()};
        const vector<GridGho>& subGrid {g.getSubG()};
        for (auto id = 1; id < mpi().getCnt(); ++id) {  // get slave buffer
            valarray<T> srcBuff(subGrid[id].N);
            mpi().Recv(&srcBuff[0], srcBuff.size() * sizeof(T), MPI_BYTE, tag,
                    id);
            glbBuff[subGs[id]] = srcBuff;
        }
        Save(glbBuff);          // Save slave buffer
    }
    else {
        mpi().Send(&lclBuff[0], lclBuff.size() * sizeof(T), MPI_BYTE, tag, 0);
    }

}

template<typename T>
void FileRaw<T>::Save(size_t glbN, valarray<T>& lclBuff, const gslice &lclGs,
                      int tag) const
{
    GsliceOpen srcGso(lclGs);          // prepare local gslice for MPI syncorize
    MPI_Datatype gsPacked;
    void* pFirst = mpi().Pack(&gsPacked, &srcGso);

    if (mpi().getId() == 0){
        valarray<T> glbBuff(glbN);            // initialize global buffer
        glbBuff[lclGs] = lclBuff;             // copy master localBuffer
        for (auto i = 1; i < mpi().getCnt(); ++i) {
            mpi().Recv(pFirst, 1, gsPacked, tag);     // gather slave gslice
            int src = mpi().getSrc();
            gslice srcGs = srcGso.getGs();

            // caculate slave length, then initialize slave buffer
            size_t length = 1;
            for (size_t a : srcGs.size()) length *= a;
            valarray<T> srcBuff(length);
            // gather slave buffer
            mpi().Recv(&srcBuff[0], srcBuff.size() * sizeof(T), MPI_BYTE, tag,
                    src);
            glbBuff[srcGs] = srcBuff;
        }
        Save(glbBuff);          // save global buffer
    }
    else {
        mpi().Send(pFirst, 1, gsPacked, tag, 0);
        mpi().Send(&lclBuff[0], lclBuff.size() * sizeof(T), MPI_BYTE, tag, 0);
    }
}

template<typename T> inline
valarray<T> FileRaw<T>::Load(size_t lenMp, const GridMc &g, int tag) const
{
    size_t lclSz = g.LG.N * lenMp;
    valarray<T> lclBuff(lclSz);

    if (mpi().getId() == 0){
        // Load all contents into global buffer
        valarray<T> glbBuff = Load(g.N * lenMp);
        size_t size = lclSz;
        size_t start = 0;
        // load data for master core
        lclBuff = glbBuff[slice(start, size, 1)];
        const vector<GridGho>& subG {g.getSubG()};
        for (auto id = 1; id < mpi().getCnt(); ++id) {
            start += size;                          // update data start pointer
            size = subG[id].N * lenMp;
            mpi().Send(&glbBuff[start], size * sizeof(T), MPI_BYTE, tag, id);
        }
    }
    else {
        mpi().Recv(&(lclBuff[0]), lclSz * sizeof(T), MPI_BYTE, tag, 0);
    }
    return lclBuff;
}

template<typename T> inline
valarray<T> FileRaw<T>::Load(size_t glbN, size_t lclN, int tag) const
{
    valarray<T> lclBuff(lclN);

    if (mpi().getId() == 0){
        // Load all contents into global buffer
        valarray<T> glbBuff = Load(glbN);
        size_t size = lclN;
        size_t start = 0;
        // load data for master core
        lclBuff = glbBuff[slice(start, size, 1)];
        for (auto id = 1; id < mpi().getCnt(); ++id) {
            start += size;                          // update data start pointer
            mpi().Recv(&size, 1, MPI_UNSIGNED_LONG, tag, id);
            mpi().Send(&glbBuff[start], size * sizeof(T), MPI_BYTE, tag, id);
        }
    }
    else {
        mpi().Send(&lclN, 1, MPI_UNSIGNED_LONG, tag, 0);
        mpi().Recv(&(lclBuff[0]), lclN * sizeof(T), MPI_BYTE, tag, 0);
    }
    return lclBuff;
}

template<typename T>
void FileRaw<T>::Save(size_t glbN, valarray<T>& lclBuff, int tag) const
{
    if (mpi().getId() == 0){
        valarray<T> glbBuff(glbN);      // create global buffer.
        size_t size = lclBuff.size();
        size_t start = 0;
        // copy master core local buffer
        glbBuff[slice(start, size, 1)] = lclBuff;
        for (auto id = 1; id < mpi().getCnt(); ++id) {
            start += size;
            mpi().Recv(&size, 1, MPI_UNSIGNED_LONG, tag, id);
            mpi().Recv(&glbBuff[start], size * sizeof(T), MPI_BYTE, tag, id);
        }
        Save(glbBuff);          // save global buffer
    }
    else {
        size_t size = lclBuff.size();
        mpi().Send(&size, 1, MPI_UNSIGNED_LONG, tag, 0);
        mpi().Send(&(lclBuff[0]), size * sizeof(T), MPI_BYTE, tag, 0);
    }
}

template<typename T>
void FileRaw<T>::Save(valarray<T> &lclBuff, int tag) const
{

    if (mpi().getId() == 0){
        // gether local size.
        va_idx_t lclSize(mpi().getCnt());
        lclSize[0] = lclBuff.size();
        for (auto id = 1; id < mpi().getCnt(); ++id) {
            size_t size = 0;
            mpi().Recv(&size, 1, MPI_UNSIGNED_LONG, tag);
            lclSize[mpi().getSrc()] = size;
        }
        size_t totSize = std::accumulate(begin(lclSize), end(lclSize),
                                         static_cast<size_t>(0));

        mpi().Barrier();

        valarray<T> glbBuff(totSize);      // create global buffer.
        // copy master local buffer
        glbBuff[slice(0, lclSize[0], 1)] = lclBuff;
        // copy slave local buffer
        size_t start = lclSize[0];
        for (auto id = 1; id < mpi().getCnt(); ++id) {
            mpi().Recv(&glbBuff[start], lclSize[id] * sizeof(T),
                       MPI_BYTE, tag, id);
            start += lclSize[id];
        }
        Save(glbBuff);          // save global buffer
    }
    else {
        size_t size = lclBuff.size();
        mpi().Send(&size, 1, MPI_UNSIGNED_LONG, tag, 0);
        mpi().Barrier();
        mpi().Send(&(lclBuff[0]), size * sizeof(T), MPI_BYTE, tag, 0);
    }
}

} // namespace IO

} // namespace DRL

#endif // FILERAW_H
