/*! \file simulator.h
 * \author Yu Zhou
 */
#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "roi.h"
#include "fileprm.h"
#include "filevti.h"
#include "lbcfg.h"
#include "iostream"
#include "gridmc.h"

namespace DRL {

namespace LBM {

using namespace GRID;
using namespace IO;
using std::string;
using std::cout;
using std::ostream;

/*!
 * \brief The Simulator class is the common interface of all Lattice Boltzmann
 *  fluid simulators.
 *
 * This class also provides shared attributes and helper functions for all
 *  derived simulators.
 *
 * Major steps to use the simulator:
 *  1. Construct the simulator with parameters.
 *  2. Run simulator.
 */
class Simulator
{
protected:
    /// \name Constants and parameters.
    ///@{
    const FilePrm& p;                           //!< Parameters.
    const bool publishedMrtSpectrum;

    enum EQUILFN{
        SRT, MRT
    };

    ///@}

    /// \name Grid and Lattice attributes.
    ///@{
    GridMc     g;      //!< The Grid.
    va_u8_t    raw;    //!< Raw core data (local).
    map_idx_t  gIdx;   //!< Grid indexed by cell value.
    va_idx_t   ibI;    //!< Inner boundary index - 0 for NOT boundary.
    veve_val_t ibNorm; //!< Inner boundary norm direction, start from 1.
    ve_val_t   ibSinCA;//!< Sin() value for contact angle.
    ve_val_t   ibCosCA;//!< Cosin() value for contact angle.
    veve_idx_t ibS;    //!< Solid direction of inner boundary lattice.
    veve_idx_t ibA;    //!< Axisymmatric mirror direction of inner boundary.
    veve_idx_t ibF[2]; //!< Half permeable film for phase 0 and 1.

    ///@}

    /// \name Controllers
    ///@{
    FileVti output;    //!< \brief Output file for saving macroscale veriables
                       //!<  \f$ \boldsymbol{u} \f$ and \f$ \rho \f$ as .vti.

    ///@}


    /// \name Model parameters
    ///@{
    va_val_t wi;   /*!<
                     * \brief The weights \f$ w_i \f$.
                     *
                     * In D1Q3 model:
                     * - 4/6 for the rest particles (i = 0);
                     * - 1/6 for i=1,2.
                     *
                     * In D2Q9 model:
                     * - 4/9 for the rest particles (i = 0);
                     * - 1/9 for i=1,2,3,4;
                     * - 1/36 for i=5,6,7,8.
                     *
                     * In D3Q15 model:
                     * - 1/8 for the rest particles (i = 0);
                     * - 1/8 for i=1,2,3,4,5,6;
                     * - 1/64 for i=7,8,9,10,11,12,13,14.
                     *
                     * In D3Q19 model:
                     * - 1/3 for the rest particles (i = 0);
                     * - 1/18 for i=1,2,3,4,5,6;
                     * - 1/36 for i=7,8,9,10,11,12,13,14,15,16,17,18.
                     */
    va_val_t wiT;  /*!<
                     * \brief Temp \f$ w_{i,temp} = 6 \frac{w_i}{c} \f$.
                     *
                     * \bug In D3Q15 model, \f$ c_s = \sqrt{\frac{3}{8}c} \f$.
                     *  (according to Zou, He), as a result, factor \f$
                     *  \frac{16}{3} \f$ should be applied instead of 6. \n
                     *  Only compabilty of D2Q9 and D3Q19 has been checked.
                     */

    ///@}

public:
    /*!
     * \brief Constructor.
     * \param prm - parameter file.
     */
    Simulator(const FilePrm& prm);

    /*!
     * \brief Run simulation (common interface of all simulator).
     * \param logProc - ostream for logging output.
     */
    virtual void Run(ostream& logProc = cout) = 0;

protected:
    // MPI_Wtime categories reported separately from construction/process time.
    // Initialization excludes initial field/checkpoint I/O; numerical work
    // includes streaming, MPI exchange, boundaries, macro variables, collision
    // and finite-state checks. Category maxima need not sum to the wall maximum.
    enum TimingCategory { INIT_TIME, NUMERICAL_TIME, FIELD_TIME, CHECKPOINT_TIME,
                          TIMING_CATEGORIES };
    void writeTiming(ostream& log, const double* seconds, double wall,
                     size_t completedSteps) const;

    ///\name Simulation steps.
    ///@{
    /*!
     * \brief Apply streaming process.
     *
     * Move the direction-specific densities \f$ f_i \f$ to the nearest
     *  neighbour lattice nodes.
     *
     * Calling MPI functions to transfer boundary values between each processor.
     *
     * \note Periodic boundary condition (loop) is automatically applied in the
     *  streaming step. To avoid periodic boundary, use bounce back, inlet, or
     *  outlet boundary.
     *
     * \param f - distribution function.
     * \param tsp - time stamp.
     */
    void streaming(vava_val_t* f, vava_val_t *fOld, size_t tsp);

    /*!
     * \brief Apply bounce back boundary condition.
     * \param f - velocity distribution.
     */
    void bounceBack(vava_val_t* f, const vava_val_t &fOld,
                    const veve_idx_t& ib);

    /*!
     * \brief axisymmetric
     * \param f - velocity distribution.
     */
    void axisymmetric(vava_val_t* f);

    /*!
     * \brief Sycnorize data at border of processors.
     *
     * \param data - data to be synced (size of data must contains ghost layer)
     * \param t - timestep.
     */
    template<typename T>
    void syncBorder(valarray<T> *data, size_t t);
    ///@}

    /// \name Simulation helper functions.
    ///@{
    /*!
     * \brief Initialize distribution functions.
     *
     * Apply given density to direction-specfic densities.
     *  \f[ f_i = w_i \rho \f]
     *
     * \param dist - distribution function.
     * \param density - macroscopic density.
     * \param ai - lattice index to be initialized.
     */
    void initDist(vava_val_t* dist, double density,
                  const va_idx_t& ai) const;

    va_idx_t getElem(size_t d, size_t first, size_t last);

    double getElemCount(const va_idx_t &val);

    double getKPar(std::size_t d);

    /*!
     * \brief Get matrix M for MRT.
     * \return valarray of matrix M.
     */
    va_val_t getMRTM();

    /*!
     * \brief Get diagonal matrix S for MRT.
     * \return valarray of diagonal elements of matrix S.
     */
    va_val_t getMRTS(double relaxationTime = 0.0);
    va_val_t forceForIndex(size_t idx) const;
    va_val_t forceSource(size_t i, const vava_val_t& velocity,
                         const va_val_t& force) const;
    bool updateConvergence(const vava_val_t& velocity, bool scalarsStable,
                           double tolerance, ostream& os,
                           const va_val_t* phaseIndicator = nullptr,
                           const vava_val_t* adjacentVelocity = nullptr,
                           const va_val_t* adjacentPhase = nullptr);
    void validatePopulations(const vava_val_t& populations) const;
    void prepareRestart(const string& path);
    void saveRestartMetadata(const string& path, size_t time, size_t phases);
    void validateRestart(const string& path, size_t time, size_t phases);
    vava_val_t previousReportVelocity;
    vava_val_t adjacentReportVelocity;
    va_val_t adjacentReportPhase;
    va_val_t previousPhaseIndicator;
    size_t stableReports = 0;
    bool hasReportVelocity = false;

    ///@}

private:
    ///\name Helper functions
    ///@{

    /*!
     * \brief Build boundary inside liquid.
     */
    void buildIb();

    /*!
     * \brief Generate weights \f$ w_a \f$
     * \return weights.
     */
    va_val_t getWi();

    /*!
     * \brief Get valarray of VD, GV, HF, VB, and PB.
     * \return valarray of space.
     */
    va_idx_t getSpace();
    ///@}
};

template<typename T>
void Simulator::syncBorder(valarray<T> *data, size_t t)
{
    static auto ibd = ROI::getIBd(g.LG);
    static auto obd = ROI::getOBd(g.LG);
    for (size_t q = 1; q < g.S.Q; ++q){
        valarray<T> sendData = (*data)[ibd[q]];
        T *pSend = &sendData[0];
        size_t size = sendData.size();
        valarray<T> recvData(size);
        T *pRecv = &recvData[0];
        mpi().SendRecv(pSend, pRecv, size * sizeof(T), MPI_BYTE, t,
                       g.AdjCi[g.S.Qop[q]], g.AdjCi[q]);
        (*data)[obd[g.S.Qop[q]]] = recvData;
    }

    //static auto ibd = ROI::getIBd(g.LG);
    //static auto obd = ROI::getOBd(g.LG);
    //for (size_t q = 1; q < g.S.Q; ++q){
    //    va_u8_t sendData = (*data)[ibd[q]];
    //    uint8_t *pSend = &sendData[0];
    //    size_t size = sendData.size();
    //    va_u8_t recvData(size);
    //    uint8_t *pRecv = &recvData[0];
    //    mpi().SendRecv(pSend, pRecv, size, MPI_UINT8_T, t,
    //                   g.AdjCi[g.S.Qop[q]], g.AdjCi[q]);
    //    (*data)[obd[g.S.Qop[q]]] = recvData;
    //}
}

} // namespace LBM

} // namespace DRL

#endif // SIMULATOR_H
