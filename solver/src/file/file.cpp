/*! \file file.cpp
 * \author Yu Zhou
 */
#include "file.h"

namespace DRL {

namespace IO {

File::File(const string &filename)
{
    Open(filename);
}

void File::Open(const string &filename)
{
    fname = filename;
}

bool File::isValid()
{
    FILE* file = fopen(fname.c_str(), "r");
    if (file != 0){
        fclose(file);
        return true;
    }
    return false;
}

} // namespace IO

} // namespace DRL
