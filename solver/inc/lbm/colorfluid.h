/*! \file colorfluid.h
 * \author Yu Zhou
 */
#ifndef COLORFLUID_H
#define COLORFLUID_H

#include "simulator.h"

namespace DRL {

namespace LBM {

/*!
 * \brief The ColorFluid class implements the color-fluid binary flow.
 *
 * Red and blue distributions functions \f$ f_i^R \f$ and \f$ f_i^B \f$ are
 *  introduced to represents two immisible fluids. Each phase undergoes the
 *  following streaming and collision scheme:
 *  \f[\underbrace{f_i^k(\boldsymbol{x} + \boldsymbol{e}_i\Delta t, t+ \Delta t)
 *    = f_i^k(\boldsymbol{x},t)}_{Streaming}
 *      + \underbrace{\Omega_i^k(\boldsymbol{x},t)}_{Collision} \f]
 *  Where \f[ \Omega_i^k = (\Omega_i^k)^3 \left[ (\Omega_i^k)^1 + (\Omega_i^k)^2
 *  \right] \f]
 *  - \f$ (\Omega_i^k)^1 \f$ is the BGK collision operator
 *  - \f$ (\Omega_i^k)^2 \f$ is a two-phase collision operator(perturbation
 *    step) which contributes to the mixed interfacial region and generates an
 *     interfacial tension.
 *  - \f$ (\Omega_i^k)^3 \f$ is recoloring step, which mimics the phase
 *    segregation and keeps the interface sharp.
 */
class ColorFluid : public Simulator
{
    friend class NumericalAudit; // Inspect actual lattice operators in the verification executable.
private:
    /// \name Constants
    ///@{
    static constexpr size_t NP = 2;  //!< \brief Number of phase = 2.
    const double SOLID_DENS_N;       //!< default densN 86 for solid.

    ///@}
    ///
    /// \name Simulation variables
    ///@{
    vava_val_t   f[NP];     /*!<
                                * \brief Distribution function / direction-
                                *  specific densities for each phase \f$
                                *  f_i^k \f$ (\f$ mu\,lu^{-3} \f$).
                                */
    vava_val_t   fOld[NP];  /*!<
                                * \brief Distribution function / direction-
                                *  specific densities of previous step.
                                */
    vava_val_t   v;         /*!<
                                * \brief Dimension-specific velocity (\f$
                                *  lu\,ts^{-1} \f$).
                                */
    va_val_t       dens[NP];  /*!<
                                * \brief Density for each phase
                                *  (\f$ mu\,lu^{-3} \f$).
                                */
    va_val_t       densN;     /*!<
                                * \brief The order parameter of binary
                                *  mixture. \f$ \rho^N = \frac{\rho^r -
                                *  \rho^b}{\rho^r + \rho^b} \f$
                                *
                                * The value of \f$ \rho^N = 1, -1, 0 \f$
                                *  correspondent to a purely red fluid, a
                                *  purely blue fluid, and the interface
                                *  respectively.
                                *
                                * \note Density ratio of one layer in
                                *  adjacent processor is needed. As a result
                                *  the size of densRto should be at least 1
                                *  layer larger than local grid size (use
                                *  lglGrid(wgl)). Addintional mpi process
                                *  should also be applied to syncorize the
                                *  adjacent densRto.
                                */
    va_val_t       densT;   /*!<
                                * \brief Sum of both phases' densitys.
                                *  \f$ \rho = \rho^r + \rho^b \f$
                                */
    va_val_t       tau;       /*!<
                                * \brief In-situ relaxation time \f$ \tau
                                *  \f$.
                                *
                                * To account for unequal viscosities of both
                                * fluids, use a linear interpolation to
                                * determin the relaxation time \f[ \tau =
                                * \frac{1+\rho^N}{2}\tau_R + \frac{1-\rho^N}
                                * {2}\tau_B \f].
                                */

    ///@}

    /// \name Helper matrixes for contact angle calculation.
    ///@{
    vava_val_t TibNorm;    //!< Rotate matrix for axisX to BD norm.
    vava_val_t TibNormInv; //!< Invert rotate matrix for axisX to BD norm.
    vava_val_t wallGradientInv;
    void calcWallGradientInverse();

    ///@}

    /// \name Model parameters
    ///@{
    va_val_t ai[NP];    /*!<
                          * \brief \f$ \phi^k_i \f$ is determined by \f$
                          *  \alpha_k \f$, which is a free parameter
                          *   determins the speed of sound.
                          *
                          * For D2Q9 model:
                          *  \f[ \phi_i^k= \left\{ \begin{array}{l@{\quad
                          *  \quad}l} \alpha_k, & i=0 \\ (1-\alpha_k)/5, &
                          *  i=1,2,3,4. \\ (1-\alpha_k)/20, & i=5,6,7,8.
                          *  \end{array} \right. \f]
                          *
                          * For D3Q19 model:
                          *  \f[ \phi_i^k= \left\{ \begin{array}{l@{\quad
                          *  \quad}l} \alpha_k, & i=0 \\ (1-\alpha_k)/12, &
                          *  i=1,2,...,6. \\ (1-\alpha_k)/24, & i=7,8,...,
                          *  18. \end{array} \right. \f]
                          *
                          * Where \f$ \alpha \f$ is a free parameter. For
                          *  D2Q9 model: \f$ (c_s^k)^2 =
                          *  \frac{3(1-\alpha_k)}{5} \f$; for D3Q19 model:
                          *  \f$ (c_s^k)^2 = \frac{1 - \alpha_k}{2} \f$.
                          *  And \f[ \lambda = \frac {\rho_R}{\rho_B} =
                          *  \frac{1-\alpha_B}{1- \alpha_R} \f]
                          *
                          * \note \f$ p_k = \rho_k (c_s^k)^2 \f$
                          */
    ///@}

public:
    /*!
     * \brief Constructor.
     * \param prm - parameter file.
     */
    ColorFluid(const FilePrm& prm);

    /*!
     * \brief Run simulation.
     *
     * The simulation is carried out in the following flow:
     * \dot
     *  digraph flow_chart{
     *    node [shape=box];
     *    start -> tStart;
     *    start   [ shape=hexagon, label="Start" ];
     *    tStart -> loadRst [ label = "N" ];
     *    tStart -> init [ label = "Y" ];
     *    tStart  [ shape=diamond, label="t=0?" ];
     *
     *    subgraph cluster0{
     *      label="Restart";
     *      loadRst -> calcM1;
     *      loadRst [ label="Load restart"];
     *      calcM1  [ label="Calculate macroscopic\n variables"];
     *    }
     *
     *    subgraph cluster1{
     *      label="Initialize";
     *      init -> calcM2 -> svOpt1;
     *      init    [ shape=box, label="Apply initial values"];
     *      calcM2  [ label="Calculate macroscopic\n variables"];
     *      svOpt1  [ label="Save output file"];
     *    }
     *
     *    { svOpt1 calcM1 } -> "t++";
     *
     *    subgraph cluster2{
     *      label="Main loop";
     *      subgraph cluster3{
     *        label="Core";
     *        style=filled;
     *        color=lightgray;
     *        node [ style=filled; color=white ];
     *        Streaming -> "Calculate macroscopic variables" -> Collision -> BC;
     *        BC [ label= "Apply boundary Conditions"];
     *      }
     *      "t++" -> Streaming;
     *      BC -> isSaveOpt;
     *      isSaveOpt -> "Save output" [ label="Y" ];
     *      isSaveOpt -> isSaveRst [ label="N" ];
     *      "Save output" -> isSaveRst;
     *      isSaveRst -> "Save restart" [ label="Y" ];
     *      isSaveRst -> isEnd [ label="N" ];
     *      "Save restart" -> isEnd;
     *      isSaveOpt [ shape=diamond; label="Save\n output?" ];
     *      isSaveRst [ shape=diamond; label="Save\n restart?" ];
     *    }
     *    isEnd -> End [ label="Y" ];
     *    isEnd -> "t++" [ label="N" ];
     *    isEnd [ shape=diamond; label="t=tEnd?" ];
     *    End [ shape=hexagon ];
     *  }
     * \enddot
     *
     * \param logProc - ostream for logging output.
     */
    virtual void Run(std::ostream &logProc = cout) override;

private:
    /// \name Major simulation steps
    ///@{
    /*!
     * \brief Apply initial values.
     *
     * Initial direction-specfic densities of both fluid are applied to void,
     *  gravity void, pressure boundary, velocity boundary, and bounceback
     *  boundary.
     *
     * \note Bounceback boundary should also have distribution functions,
     *  otherwise density in void will occupy it after streaming, the total
     *  mass in void will be reduced.
     */
    void initialize();

    /*!
     * \brief Calculate macroscopic variables.
     *
     * The macroscopic fluid density: \f[
     * \rho = \sum\limits_{a=0}^{Q} f_i \f]
     * The macroscopic velocity \f$ \boldsymbol{u} \f$ is an average of the
     *  microscopic velocities \f$ \boldsymbol{e}_i \f$ weighted by the
     *  directional densities \f$ f_i \f$: \f[
     *  \boldsymbol{u} = \frac{1}{\rho} \sum\limits_{a=0}^Q f_i \boldsymbol{e}_i
     *  \f]
     *
     * Applis to grids indexed "PB", "VB", "BB", and "VD".
     *
     * \param t - timestep, mpi tag for this operation.
     */
    void calcMacroVar(size_t t);

    /*!
     * \brief The collision step.
     *
     * Calculate \f$ f^{eq} \f$ in tems of \f$ u \f$, then perform the key step
     *  as described in the SCSP class description.
     *
     * The collosion step only applied to grid indexed "VD", "GV", and "HF".
     *
     * \note For "GV", the gravity is incorporated in the calculation of \f$
     *  f^{eq} \f$, where \f$ \boldsymbol{u}^{eq} \f$ used in computation of \f$
     *  f^{eq} \f$ is introduced as: \f[
     *  \boldsymbol{u}^{eq} = \boldsymbol{u} + \Delta \boldsymbol{u}
     *                  = \boldsymbol{u} + \frac{\tau \boldsymbol{F}}{\rho} \f]
     *  where \f[
     *  \boldsymbol{F} = m \boldsymbol{a}
     *                 = m \frac{\mathrm{d}\boldsymbol{u}}{\mathrm{d}t}
     *                                                               \Rightarrow
     *  \Delta \boldsymbol{u} = \frac{\tau \boldsymbol{f}}{\rho} \f]
     *
     * \param idx - index of the grid to be applied to.
     */
    void collision();

    /*!
     * \brief Apply Zou & He boundary conditions.
     *
     * The following boundary conditions are applied:
     * 1. Bounceback Boundary condition (for solid-fluid contact face)
     * ------------------------------
     *  - Bounceback for grides indexed "BB". Using "mid-plane" bounceback
     *    scheme in which the directional density is stored inside solids,
     *    bounce back, and re-emerge at the next time step. In programming: \f$
     *    switch(f_i, opposit(f_i)) \f$.
     *
     * 2. Zou and He boundary condition (for pressure and flux)
     * ----------------------------------------------------
     * ### 2.1 Priciples###
     *  - Von Neumann (flux) boundaries for grids indexed "VB", given velocity
     *    \f$ \boldsymbol{u}_0 \f$. The following equations should be satisfied
     *    for monentum and mass conversation: \f[
     *    \rho = \sum\limits_{a}f_a, \qquad
     *    \boldsymbol{u}_0 = \frac{1}{\rho}\sum\limits_{a}f_a\boldsymbol{e}_a\f]
     *    And by applying Zou and He boundary condition for 2D: \f[
     *    \begin{eqnarray}
     *    f_i - f_i^{eq} = f_{op(i)} - f_{op(i)}^{eq} &\quad
     *                     (\underbrace{i = 1, 2, 3, 4}_{for~2D})&\quad
     *                     (\underbrace{i = 1, 2, \ldots, 6}_{for~3D}) \\
     *    f_i - f_i^{eq} = f_{op(i)} - f_{op(i)}^{eq} +
     *                     \frac{1}{2}\boldsymbol{e}_i\boldsymbol{\delta} &\quad
     *                     (\underbrace{i = 5, 6, 7, 8}_{for~2D})&\quad
     *                     (\underbrace{i = 7, 8, \ldots, 18}_{fir~3D})
     *    \end{eqnarray}\f]
     *    \f$ op(i) \f$ for opposit direction. Note that \f$ \boldsymbol{\delta}
     *    \f$ is added to correct directional monenta respectively, but leave
     *    the dirction normal to the boundary _unchanged_.
     *    To make it clear, take the D2Q9 scheme defined in the GridQ class for
     *    example , say, west boundary (fluid direction x+), \f$ f_1,\;f_5,\;
     *    f_8,\;\rho \f$ are unknowns. The following equations can be derived:
     *    \f[ \begin{eqnarray}
     *    f_1 &=& f_3 + (f_1^{eq} - f_3^{eq})\\
     *        &=& f_3 + \frac{2}{3}\rho \frac{u_x}{c^2} &\quad(1)\\
     *    f_5 &=& f_7 + (f_5^{eq} - f_7^{eq}) + \frac{1}{2}e_{5y}\delta_y\\
     *        &=& f_7 + \frac{1}{6}\rho \frac{u_x + u_y}{c^2}
     *                + \frac{1}{2}(\frac{2}{3}\rho u_y - f_2 + f_4) &\quad(2)\\
     *    f_8 &=& f_6 + (f_8^{eq} - f_6^{eq}) + \frac{1}{2}e_{8y}\delta_y\\
     *        &=& f_6 + \frac{1}{6}\rho \frac{u_x - u_y}{c^2}
     *                - \frac{1}{2}(\frac{2}{3}\rho u_y - f_2 + f_4) &\quad(3)\\
     *    \rho &=& f_0 + f_1 + f_2 + f_3 + f_4 + f_5 + f_6 + f_7 + f_8
     *             &\quad(4)\\
     *    \rho u_x &=& f_1 + f_5 + f_8 - f_3 -f_7 - f_6 &\quad(5)\\
     *    \rho u_y &=& f_2 + f_5 + f_6 - f_4 - f_7 - f_8 &\quad(6)
     *    \end{eqnarray}  \f]
     *    Note that equition \f$ (5) = (1) + (2) + (3) \f$, and equition \f$ (6)
     *    = (2) - (3) \f$. Four independent equitions for four unknown variables
     *    to solve. For D3Q19 scheme, say, west boundary (fluid direction x+),
     *    \f$ f_1,\;f_7,\;f_8,\;f_{15},\;f_{16},\;\rho \f$ are unknowns. The
     *    following equations can be derived:
     *    \f[ \begin{eqnarray}
     *    f_1 &=& f_2 + (f_1^{eq} - f_2^{eq})\\
     *        &=& f_2 + \frac{1}{3}\rho \frac{u_x}{c^2} &\quad(7)\\
     *    f_7 &=& f_9 + (f_7^{eq} - f_9^{eq})
     *                + \frac{1}{2}e_{7y}\delta_y \\
     *        &=& f_9 + \frac{1}{6}\rho \frac{u_x + u_y}{c^2}
     *                + \frac{1}{2}(\frac{2}{3}\rho u_y - f_3 + f_4 - f_{11}
     *                              + f_{13} - f_{12} + f_{14}) &\quad(8) \\
     *    f_8 &=& f_{10} + (f_8^{eq} - f_{10}^{eq})
     *                   + \frac{1}{2}e_{8y}\delta_y \\
     *        &=& f_{10} + \frac{1}{6}\rho \frac{u_x - u_y}{c^2}
     *                   - \frac{1}{2}(\frac{2}{3}\rho u_y - f_3 + f_4 - f_{11}
     *                                 + f_{13} - f_{12} + f_{14}) &\quad(9) \\
     *    f_{15} &=& f_{17} + (f_{15}^{eq} - f_{17}^{eq})
     *                      + \frac{1}{2}e_{15z}\delta_z \\
     *           &=& f_{17} + \frac{1}{6}\rho \frac{u_x + u_z}{c^2}
     *                      + \frac{1}{2}(\frac{2}{3}\rho u_z - f_5 + f_6
     *                                    - f_{11} + f_{13} + f_{12} - f_{14})
     *                                                            &\quad(10) \\
     *    f_{16} &=& f_{18} + (f_{16}^{eq} - f_{18}^{eq})
     *                      + \frac{1}{2}e_{16z}\delta_z \\
     *           &=& f_{18} + \frac{1}{6}\rho \frac{u_x - u_z}{c^2}
     *                      - \frac{1}{2}(\frac{2}{3}\rho u_z - f_5 + f_6
     *                                    - f_{11} + f_{13} + f_{12} - f_{14})
     *                                                            &\quad(11) \\
     *    \rho &=& f_0 + f_1 + f_2 + f_3 + f_4 + f_5 + f_6 + f_7 + f_8 + f_9 +
     *             f_{10} + f_{11} + f_{12} + f_{13} + f_{14} + f_{15} +
     *             f_{16} + f_{17} + f_{18} &\quad(12)\\
     *    \rho u_x &=& f_1 + f_7 + f_8 + f_{15} + f_{16} -
     *                 f_2 - f_9 - f_{10} - f_{17} - f_{18} &\quad(13)\\
     *    \rho u_y &=& f_3 + f_7 + f_{10} + f_{11} + f_{12} -
     *                 f_4 - f_8 - f_9 - f_{13} - f_{14} &\quad(14)\\
     *    \rho u_z &=& f_5 + f_{11} + f_{14} + f_{15} + f_{18} -
     *                 f_6 - f_{12} - f_{13} - f_{16} - f_{17} &\quad(15)\\
     *    \end{eqnarray}  \f]
     *    Note that equition \f$ (13) = (7) + (8) + (9) + (10) + (11) \f$,
     *    equition \f$ (14) = (8) - (9) \f$, and equition \f$ (15) = (10) - (11)
     *    \f$. Six independent equitions for six unknown variables to solve.
     *  - Dirichlet (pressure) boundaries for grides indexed "PB". Note that
     *    speciifying pressure is equivalent to specifying density since there
     *    is an equation of state relating them directly. The relationship is
     *    simply \f$ P = \rho RT = \rho c_s^2 \f$ The derivation of pressure
     *    boundary condition is similar to flux boundary. Unknown variables are
     *    computed from equation (1 - 15).
     *
     * ###2.2 Unified boundary drived from above###
     * For simplicity, suppose flow direction is always parallel to x or y or z
     * axis, which means the velocity parallel to boundary ( normal to the flow
     * direction, in other words) is 0. For example flow from \f$ yz \f$ plane
     * has \f$ u_y = 0, ~u_z = 0 \f$. In this case, we only need to determin
     * \f$ u_x \f$ and unknown \f$ f_i \f$.
     *  - Calculation of unknown velocity and density:
     *   \f[ \rho = \frac
     *   {2\sum{\boldsymbol{f_{opp}}} + \sum{\boldsymbol{f_{vertical}}}}
     *   {1-\boldsymbol{u}} \f]
     *  - Calculation of unknown distribution function:
     *   \f[ f_{i} = f_{opp} + 6w_i\rho\frac{\boldsymbol{u}}{c} +
     *   \frac{1}{2}\sum{f_{vertical}e_i}\f]
     *   \f$ f_{opp} \f$ - opposit direction \n
     *   \f$ f_{vertical} \f$ - directions vertical to the flow direction.
     *
     * \warning This function may produce wrong results by having the boundary
     *  at the border lattice in a local grid.
     */
    void zouHeBoundary();

    /*!
     * \brief Apply body force.
     *
     * \see collision()
     */
    void bodyForce();

    /*!
     * \brief Apply half permeable film.
     *
     * Half permeable film only let allowed phase pass, and apply bounceback to
     *  the other phase.
     */
    void clearHalfFilm();
    ///@}

    ///\name I/O operations
    ///@{
    /*!
     * \brief Load restart file.
     *
     * Get the seralized data and recover it to the local distribution function.
     * \sa saveRst() for contrary operation.
     * \param t - timestep, mpi tag for this operation.
     */
    void loadRst(size_t t);

    /*!
     * \brief Save restart file.
     *
     * Seralize local distribution function and save it.
     * \sa loadRst() for contrary operation.
     * \param t - timestep, mpi tag for this operation.
     */
    void saveRst(size_t t);

    /*!
     * \brief Report simulation status.
     * \param os - output streaming.
     * \param r0 - density of phase 0.
     * \param r1 - density of phase 1.
     * \param v - velocity.
     * \return true if satisfied ending condition.
     */
    bool report(ostream &os, const va_val_t& r0, const va_val_t& r1,
                const vava_val_t &v);
    ///@}

    /// \name Helper functions
    ///@{
    /*!
     * \brief Calculate translation (rotate) matrix for inner boundary normal.
     *
     * X axis multiplied by this matrix is rotated to the normal of inner
     *  boundary.
     *
     * Also calculate the inverse translation matrix.
     */
    void calcTibNorm();

    /*!
     * \brief Initialize \f$ \Phi_i \f$.
     */
    void initPhi();

    /*!
     * \brief Generate weights \f$ b_a \f$
     * \return array of weights.
     */
    va_val_t getBi();
    ///@}
};

} // namespace LBM

} // namespace DRL

#endif // COLORFLUID_H
