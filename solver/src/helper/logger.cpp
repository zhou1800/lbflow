/*! \file logger.cpp
 * \author Yu Zhou (zhouyu@live.com)
 */
#include "logger.h"
#include <sstream>
#include <iomanip>

namespace DRL {

namespace HLP {

using std::ostringstream;
using std::setfill;
using std::setw;
using std::setprecision;
using std::setiosflags;

Logger::Logger(std::ostream &os, bool master)
    : out(os)
    , timer(mpi().getTime()), isMaster(master), isAutoRst(false)
{

}

void Logger::Log(const string &s, bool showTimer)
{
    if (isMaster && mpi().getId() != 0)
        return;
    out << s;
    if (showTimer)
        out << "Time elapsed: " << hTime(mpi().getTime() - timer) << std::endl;
    if (isAutoRst)
        timer = mpi().getTime();
    out.flush();
}

void Logger::Reset(bool autoReset)
{
    isAutoRst = autoReset;
    timer = mpi().getTime();
}

void Logger::Pause()
{
    if (isPaused) return;
    isPaused = true;
    paused = mpi().getTime();
}

void Logger::Continue()
{
    if (!isPaused) return;
    isPaused = false;
    timer += mpi().getTime() - paused;
}

std::string Logger::hTime(double time)
{
    const size_t sInMinute = 60;
    const size_t sInHour = 3600;
    const size_t sInDays = 86400;
    size_t iTime = static_cast<size_t>(time);
    size_t days = iTime / sInDays;
    size_t hh = (iTime % sInDays) / sInHour;
    size_t mm = (iTime % sInHour) / sInMinute;
    size_t ss = iTime % sInMinute;
    double sLeft = time - static_cast<double>(iTime) + ss;
    ostringstream os;
    os.precision(6);
    os << days << "d,"
        << setw(2) << setfill('0') << hh << ":"
        << setw(2) << setfill('0') << mm << ":"
        << setw(9) << setfill('0') << setiosflags(std::ios::fixed) << sLeft;
    return os.str();
}

} // namespace HLP

} // namespace DRL

