/*! \file histogram.h
 * \author Yu Zhou (zhouyu@live.com)
 */

#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include "imgproc.h"

namespace DRL {

namespace PIPT {

/*!
 * \brief The Histogram class calculate histogram of given raw file.
 *
 * Firstly initialize the multi-core grid, then load the raw file, and calculate
 *  the histogram. Secondly, calculate the accumulates and percentage of the
 *  histogram. Finally save the results as a csv file.
 */
class Histogram : public ImgProc
{
private:
    /// \name Histogram and statistics.
    ///@{
    va_u64_t hist;     /*!<
                         * \brief Histogram.
                         *
                         * Initialized in calHist(). Contains global histogram
                         *  on master processor, and local on slave processors.
                         */
    va_val_t accu;     ///< Accumulates - Initailized in calSum().
    va_val_t pct;      ///< Percentage - Initialized in calSum().

    ///@}

public:
    /// \name Constructors and destructor.
    ///@{
    /*!
     * \brief Construct with existing GridMc.
     * \param g - pre-initialized global grid.
     */
    Histogram(const GridMc& g);

    /*!
     * \brief Construct and create new GridMc
     * \param sz - dimension size of the grid.
     * \param q - number of directions.
     * \param ghoSz - ghost layer size (default = 1).
     */
    Histogram(const va_idx_t& sz, size_t q, size_t ghoSz = 1);
    ///@}

    /// \name Main operation.
    ///@{
    /*!
     * \brief Calculate histogram.
     *
     * Calculate histogram, accumulates and percentage.
     *
     * \sa calStatistics()
     * \param buff - raw buffer to be calculated.
     * \return histogram.
     */
    template <typename T>
    const va_u64_t& CalHist(const valarray<T>& buff);

    /*!
     * \brief Save histogram file.
     * \param logPath - file path to be saved.
     * \param logFname - file name to be saved.
     * \return actual saved file path name.
     */
    string SaveHistFile(const string &logPath, const string &logFname);
    ///@}

    /// \name Get functions.
    ///@{
    /*!
     * \brief Get histogram.
     * \return histogram.
     */
    const va_u64_t& getHist() const { return hist; }

    /*!
     * \brief Get accumulates.
     * \return accumulates.
     */
    const va_val_t& getSumAccu() const { return accu; }

    /*!
     * \brief Get percentage.
     * \return percentages.
     */
    const va_val_t& getSumPct() const { return pct; }

    double getPorosity(size_t threshold) const
    { return pct[threshold >= pct.size() ? pct.size() - 1 : threshold]; }
    ///@}

private:
    /// \name Helper functions.
    ///@{
    /*!
     * \brief Calculate accumulation and percentage.
     */
    void calStatistics();

    ///@}
};

template <typename T> inline
const va_u64_t& Histogram::CalHist(const valarray<T> &buff)
{
    // Calculate local histogram
    hist = va_u64_t((uint64_t)0, (size_t)(std::pow(256, sizeof(T))));
    for (T a : buff)
        hist[a]++;

    // Accomulate local histogram and distribute to all processors.
    va_u64_t sum(hist.size());
    mpi().AllReduce(&hist[0], &sum[0], sum.size(), MPI_UINT64_T, MPI_SUM);
    hist = sum;

    // Calculate hisogram statics
    calStatistics();

    return hist;
}

} // namespace PIPT

} // namespace DRL

#endif // HISTOGRAM_H
