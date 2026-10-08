/*! \file logger.h
 * \author Yu Zhou (zhouyu@live.com)
 */
#ifndef LOGGER_H
#define LOGGER_H

#include "thread.h"
#include <string>

namespace DRL {

namespace HLP {

using std::ostream;
using std::string;

/*!
 * \brief The Logger class manage the log file.
 */
class Logger
{
private:
    ostream& out;           //!< Out stream for log.
    double timer;           //!< Timer.
    bool isMaster;          //!< Only master uses the logger.
    bool isAutoRst;         //!< Timer auto reset after each log.

    bool isPaused;          //!< Flag of timer pause.
    double paused;          //!< Paused time.

public:
    /*!
     * \brief Constructor.
     * \param os - out put stream.
     * \param master - only log from master core if true.
     */
    Logger(ostream& os, bool master = true);

    /// \name Log operation.
    ///@{
    /*!
     * \brief Log to outstream.
     *
     * Also reset timer if isAutoRst is true.
     *
     * \sa Reset() for set isAutoRst.
     * \param s - string to be logged.
     * \param showTimer - append timer to log if true.
     */
    void Log(const string& s, bool showTimer = true);
    ///@}

    /// \name Timer operation
    ///@{
    /*!
     * \brief Reset timer, also toggle auto reset flag.
     * \param autoReset - auto reset timer flag.
     */
    void Reset(bool autoReset = false);

    /*!
     * \brief Pause timer
     */
    void Pause();

    /*!
     * \brief Continue paused timer.
     */
    void Continue();
    ///@}

private:
    /*!
     * \brief Convert time (in seconds) to human readable time (HH:MM:SS.xxx)
     * \param time - time to be convert.
     * \return human readable time.
     */
    string hTime(double time);
};

} // namespace HLP

} // namespace DRL

#endif // LOGGER_H
