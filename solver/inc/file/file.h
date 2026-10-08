/*! \file file.h
 * \author Yu Zhou (zhouyu@live.com)
 */
#ifndef FILE_H
#define FILE_H

#include "defs.h"

namespace DRL {

namespace IO {

/*!
 * \brief The File class provides file manipulation interface.
 *
 * Stores a filename. Provides methods to check if the file exists.
 */
class File
{
protected:
    string fname;   ///< File name.

public:
    /*!
     * \brief Constructor.
     * \sa Open()
     * \param filename - file for input and output.
     */
    File(const string& filename);

    /// \name File states functions.
    ///@{
    /*!
     * \brief Open a file for input and output.
     * \param filename - file name.
     */
    virtual void Open(const string& filename);

    /*!
     * \brief Test if the opened file is valid.
     *
     * Try to open the file to test its existece.
     *
     * \return true if file exists.
     */
    bool isValid();
    ///@}
};

} // namespace IO

} // namespace DRL

#endif // IO_H
