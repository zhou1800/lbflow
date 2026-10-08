/*! \file filevti.h
 * \author Yu Zhou
 */
#ifndef FILEVTI_H
#define FILEVTI_H

#include "file.h"
#include <list>

namespace DRL {

namespace IO {

using std::list;

/*!
 * \brief The FileVti class save .vti format file to be open with ParaView.
 *
 * This class provides methods to write float point variables into a
 *  VTK(ParavView) xml style .vti file.
 *
 * Useage:
 *  First add variables by calling AddScalar() and AddVector(). Then call Init()
 *  to generate prefix and appendix in output data. Finally call Save() to save
 *  .vti file.
 *
 * The .vti file format:
 *
 ~~~~~~~~~~{.xml}
 * <VTKFile type="ImageData" version="1.0" byte_order="LittleEndian" header_type="UInt64">
 *   <ImageData WholeExtent="0 nx 0 ny 0 nz" Origin="0 0 0" Spacing="1 1 1">
 *     <Piece Extent="x0 x1 y0 y1 z0 z1">
 *       <PointData Scalars="s1 s2" Vectors="v1">
 *         <DataArray type="Float64" Name="s1" format="appended" offset="appOff0"/>
 *         <DataArray type="Float64" Name="s2" format="appended" offset="appOff1"/>
 *         <DataArray type="Float64" Name="v1" NumberOfComponents="3" format="appended" offset="appOff2"/>
 *       </PointData>
 *       <CellData></CellData>
 *     </Piece>
 *   </ImageData>
 *   <AppendedData encoding=\"raw\">
 *     _DATA
 *   </AppendedData>
 * </VTKFile>
 ~~~~~~~~~~
 *
 * Appended with data. A data length with type of uint64_t is added before each
 *  section of _DATA.
 *
 * Embedding in C++:
 *
 * | Line  | Name      | Distribution | Parameters        | Value         |
 * |-------|-----------|--------------|-------------------|-------------- |
 * | 1~2   | headStart | master       | nx ny nz          | glbGrid       |
 * | 3~10  | piece     | all          | x0 x1 y0 y1 z0 z1 | glbGrid       |
 * |       | piece     | all          | s1 s2 v1          | addVariable() |
 * | 11~12 | headEnd   | master       | none              | N/A           |
 * | 13    |           | all          | DATA              | addVariable() |
 * | 14~15 | foot      | master       | none              | N/A           |
 *
 */
class FileVti : public File
{
private:
    /// \name .vti file perfix and appendix
    ///@{
    string headStart;       //!< Start of .vti file head.
    string piece;           //!< Pieces in .vti file head.
    string headEnd;         //!< End of .vti file head.
    string foot;            //!< Foot of .vti file.

    ///@}

    /// \name Pointer to data to be saved
    ///@{
    map<string, va_val_t*>       scalars;   //!< Scalars to be saved.
    map<string, list<va_val_t*>> vectors;   //!< Vectors to be saved.

    ///@}
public:
    /*!
     * \brief Construct with a filename.
     * \param filename - output .vti file.
     */
    FileVti(const string &filename);

    /// \name Preparation
    ///@{
    /*!
     * \brief Add a scalar to the .vti file.
     *
     * \note Scalar name must be unique, otherwise it will override the previous one.
     *
     * \param name - scalar name.
     * \param value - pointer to scalar value.
     */
    void AddScalar(const string& name, va_val_t *value);

    /*!
     * \brief Add a vector to the .vti file.
     *
     * New added vector value will be append to the existing vectors with the
     *  same name.
     *
     * \param name - vector name.
     * \param value - pointer to the value of vector.
     */
    void AddVector(const string& name, va_val_t *value);

    /*!
     * \brief Construct the .vti file header with grid information.
     *
     * Generate prefix and appendix in .vti file.
     *
     * \param D - number of dimensions.
     * \param glbN - number of elements in global grid.
     * \param glbSz - global grid size.
     * \param glbStride - stride of global grid.
     * \param lclSz - local grid size.
     * \param lclStart - start point of local grid.
     */
    void Init(size_t D, size_t glbN,
              const va_idx_t &glbSz, const va_idx_t &glbStride,
              const va_idx_t &lclSz, const va_idx_t &lclStart);
    ///@}

    /// \name Save .vti file
    ///@{
    /*!
     * \brief Save the .vti file.
     *
     * Perform the save operation. The .vti file name is generated automatically
     *  as file name + file sequence.
     * \param seq - file sequence.
     */
    void Save(size_t seq);
    ///@}

    /// \name Clean up and reverse.
    ///@{
    /*!
     * \brief Remove all added variables and generated .vti file information.
     */
    void clean();
    ///@}
};

} // namespace IO

} // namespace DRL

#endif // FILEVTI_H
