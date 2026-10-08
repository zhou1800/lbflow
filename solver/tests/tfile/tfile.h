#ifndef TFILE_H
#define TFILE_H

#include "debugger.h"
#include "file.h"
#include "fileraw.h"
#include "filevti.h"
#include "fileprm.h"
#include <algorithm>

using namespace DRL;
using namespace DRL::TEST;
using namespace DRL::IO;
using namespace DRL::GRID;
using namespace std;

void tFile(ostream& os, const string& fname);
void tFilePrm(ostream& os, const string& fname);

template <typename T>
void tFileRaw(ostream& os, const FilePrm& prm);
template <typename T>
bool compFileRaw(const FileRaw<T>& iraw, const FileRaw<T>& oraw, size_t N);

void tFileVti(ostream& os, const FilePrm& prm);

template <typename T>
void tFileRaw(ostream& os, const FilePrm &prm)
{
    FileRaw<T> iraw(prm.getGFile());
    FileRaw<T> oraw(prm.getOutPath());
    GridMc g(prm.getD(), prm.getQ());
    if (!compFileRaw<T>(iraw, oraw, g.N))
        os << "Load1 or Save1 error." << endl;

    valarray<T> buff2 = iraw.Load(g);                 // distribute load. Load2
    valarray<T> buff3 = iraw.Load(g.N, g.LG.N, g.LGs);// distribute load. Load3
    if (!equal(begin(buff2), end(buff2), begin(buff3)))
        os << "Load2 or Load3 mismatch." << endl;

    oraw.Save(g, buff2);                              // distribute save. Save2
    if (!compFileRaw<T>(iraw, oraw, g.N))
        os << "Save2 or Load2/3 error." << endl;

    oraw.Save(g.N, buff2, g.LGs);                     // distribute save. Save3
    if (!compFileRaw<T>(iraw, oraw, g.N))
        os << "Save3 error." << endl;

    valarray<T> buff4 = iraw.Load(g, 1, MPI_ANY_TAG); // distribute load. Load4
    valarray<T> buffT = buff4[g.LG.getI2fGs()];
    oraw.Save(g, buffT);
    if (!compFileRaw<T>(iraw, oraw, g.N))
        os << "Load4 error." << endl;

    valarray<T> buff5 = iraw.Load(1, g);              // Seralized Load. Load5
    oraw.Save(g.N, buff5);                            // Seralized Save. Save4
    if (!compFileRaw<T>(iraw, oraw, g.N))
        os << "Load5 or Save4 error." << endl;

    valarray<T> buff6 = iraw.Load(g.N, g.LG.N);              // Seralized Load. Load6
    oraw.Save(buff6, MPI_ANY_TAG);                    // Seralized Save. Save5
    if (!compFileRaw<T>(iraw, oraw, g.N))
        os << "Load6 or Save5 error." << endl;

    os << "Test FileRaw Passed!" << endl;
}

template <typename T>
bool compFileRaw(const FileRaw<T> &iraw, const FileRaw<T> &oraw, size_t N)
{
    if (mpi().getId() == 0) {
        valarray<T> buff0 = iraw.Load(N);    // load origin. Load 1
        oraw.Save(buff0);                      // Save origin. Save 1
        valarray<T> buff1 = oraw.Load(N);      // load duplicate. Load 1
        if (!equal(begin(buff1), end(buff1), begin(buff0))) // compare.
            return false;   // note: only head core return false.
    }
    else
        return true;
}

#endif // TFILE_H
