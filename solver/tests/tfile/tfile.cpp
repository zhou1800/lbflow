#include "tfile.h"

void tFile(ostream &os, const string &fname)
{
    File f(fname);
    if (f.isValid())
        os << "The file " << fname << " is valid." << endl;
    else
        os << "The file " << fname << " is not valid." << endl;
}
