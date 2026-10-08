/*! \file pipt/histogram.cpp
 * \author Yu Zhou (zhouyu@live.com)
 */

#include "histogram.h"
#include <fstream>
#include <numeric>

namespace DRL {

namespace PIPT {

using namespace std;

Histogram::Histogram(const GridMc &g)
    : ImgProc(g)
    , hist(0), accu(0), pct(0)
{

}

Histogram::Histogram(const va_idx_t &sz, size_t q, size_t ghoSz)
    : ImgProc(sz, q, ghoSz)
    , hist(0), accu(0), pct(0)
{

}

string Histogram::SaveHistFile(const string& logPath, const string& logFname)
{
    if (mpi().getId() == 0) {
        char name[255];
        sprintf(name, "%s/%s.csv", logPath.c_str(), logFname.c_str());
        ofstream file(name);
        for (size_t i = 0; i < hist.size(); ++i)
            file << i << "," << hist[i] << "," << accu[i] << ","
                 << pct[i] << "," << endl;
        file.close();
        return name;
    }
    return string();
}

void Histogram::calStatistics()
{
    double total = (double)G.N / 100.0;
    accu.resize(hist.size());
    partial_sum(begin(hist), end(hist), begin(accu));
    pct = accu / total;
}

} // namespace PIPT

} // namespace DRL
