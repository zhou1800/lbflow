/*! \file filevti.cpp
 * \author Yu Zhou
 */
#include "filevti.h"
#include "fileraw.h"
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cstring>

namespace DRL {

namespace IO {

using std::runtime_error;
using std::stringstream;
using std::ostringstream;
using std::setfill;
using std::setw;

FileVti::FileVti(const string &filename)
    : File(filename)
    , headStart(), piece(), headEnd(), foot()
    , scalars(), vectors()
{
}

void FileVti::AddScalar(const std::string &name, va_val_t* value)
{
    scalars[name] = value;
}

void FileVti::AddVector(const std::string &name, va_val_t* value)
{
    vectors[name].push_back(value);
}

void FileVti::Init(size_t D, size_t glbN,
                   const va_idx_t& glbSz, const va_idx_t& glbStride,
                   const va_idx_t& lclSz, const va_idx_t& lclStart)
{
    /* write headStart */
    ostringstream os;
    if (mpi().getId() == 0) {
        os << "<VTKFile type=\"ImageData\" version=\"1.0\""
              " byte_order=\"LittleEndian\" header_type=\"UInt64\">"
           << "<ImageData WholeExtent=\" ";
        for (size_t d = D; d > 0; --d)
            os << "0 " << glbSz[d - 1] - 1 << " ";
        for (size_t d = 3 - D; d > 0; --d)
            os << "0 0 ";
        os << "\" Origin=\"0 0 0\" Spacing=\"1 1 1\">";
        headStart = os.str();
        os.str("");
    }

    /* write piece */
    os << "<Piece Extent=\"";
    for (size_t d = D; d > 0; --d)
        os << lclStart[d - 1] << " "
           << lclStart[d - 1] + lclSz[d - 1] - 1 << " ";
    for (size_t d = 3 - D; d > 0; --d)
        os << "0 0 ";
    os << "\"><PointData Scalars=\"";
    for (const auto& a : scalars)
        os << a.first << " ";
    os << "\" Vectors=\"";
    for (const auto& a : vectors)
        os << a.first << " ";
    os << "\" >";

    size_t offRaw = 0;
    size_t lclN = 1;
    for (size_t d = 0; d < D; ++d) {
        offRaw += lclStart[d] * glbStride[d] * lclN;
        lclN *= lclSz[d];
    }

    // for scalars
    size_t offset = offRaw * sizeof(double) + mpi().getId() * sizeof(uint64_t);
    for (const auto& a : scalars) {
        os << "<DataArray type=\"Float64\" Name=\"" << a.first
           << "\" format=\"appended\" offset=\"" << offset << "\"/>";
        offset += glbN * sizeof(double)
                + mpi().getCnt() * sizeof(uint64_t);
    }

    // for vectors
    if (vectors.size() > 0)
        offset = offRaw * sizeof(double) * vectors.begin()->second.size()
                + mpi().getId() * sizeof(uint64_t)
                + glbN * scalars.size() * sizeof(double)
                + mpi().getCnt() * scalars.size() * sizeof(uint64_t);
    for (const auto& a : vectors) {
        os << "<DataArray type=\"Float64\" Name=\"" << a.first
           << "\" NumberOfComponents=\"" << a.second.size()
           << "\" format=\"appended\" offset=\"" << offset << "\"/>";
        offset += glbN * sizeof(double) * 3
                + mpi().getCnt() * sizeof(uint64_t);
    }

    os << "</PointData><CellData></CellData></Piece>";
    piece = os.str();

    /* write cross and foot */
    if (mpi().getId() == 0) {
        headEnd = "</ImageData><AppendedData encoding=\"raw\">_";
        foot = "</AppendedData></VTKFile>";
    }
}

void FileVti::Save(std::size_t seq)
{
    if (scalars.size() == 0 && vectors.size() == 0)
        throw runtime_error("FileVti::Save - nothing to save.");
    if (piece.empty())
        throw runtime_error("FileVti::Save - vti file needs to be initialized"
                            " before save.");
    if (sizeof(uint64_t) != sizeof(double))
        throw runtime_error("Lattice::SaveOpt() - cannot perform conversion");

    // prepare file
    ostringstream os;
    os << fname << "." << setw(10) << setfill('0') << seq << ".vti";
    FileRaw<char> charFile(os.str(), true);
    FileRaw<double> rawFile(os.str(), true);

    // Save headStart
    if (mpi().getId() == 0)
        charFile.Save(valarray<char>(headStart.c_str(), headStart.size()));

    // Save piece
    valarray<char> piece_arr(piece.c_str(), piece.size());
    //mpi().Barrier();
    charFile.Save(piece_arr, seq);

    // Save headEnd
    if (mpi().getId() == 0)
        charFile.Save(valarray<char>(headEnd.c_str(),headEnd.size()));

    // Save Scalars
    for (const auto& a : scalars) {
        uint64_t len = a.second->size();
        va_val_t buff(len + 1);
        len *= sizeof(double);
        memcpy(&buff[0], &len, sizeof(double));
        buff[slice(1, a.second->size(), 1)] = *a.second;
        //mpi().Barrier();
        rawFile.Save(buff, seq);
    }

    // Save vectors
    for (auto& a : vectors) {
        if (a.second.size() == 0) break;
        uint64_t len = a.second.size() * a.second.front()->size();
        va_val_t buff(len + 1);
        len *= sizeof(double);
        memcpy(&buff[0], &len, sizeof(double));
        size_t start = 1;
        for (auto& b : a.second) {
            buff[slice(start, b->size(), a.second.size())] = *b;
            start++;
        }
        //mpi().Barrier();
        rawFile.Save(buff, seq);
    }

    // Save foot
    if (mpi().getId() == 0)
        charFile.Save(valarray<char>(foot.c_str(), foot.size()));
}

void FileVti::clean()
{
    headStart.clear();
    piece.clear();
    headEnd.clear();
    foot.clear();
    scalars.clear();
    vectors.clear();
}

} // namespace IO

} // namespace DRL
