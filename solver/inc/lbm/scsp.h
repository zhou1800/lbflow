/*! \file scsp.h
 * \author Yu Zhou (zhouyu@live.com)
 */
#ifndef SCSP_H
#define SCSP_H

#include "simulator.h"

namespace DRL {

namespace LBM {

/*!
 * \brief The SCSP class implements single phase single components fluid model.
 */
class SCSP : public Simulator
{
private:
    /// \name Simulation variables
    ///@{
    vava_val_t   f;        /*!<
                               * \brief Distribution function / direction-
                               *  specific densities \f$ f_i \f$ (\f$
                               *  mu\,lu^{-3} \f$).
                               */
    vava_val_t   fOld;     /*!<
                               * \brief Distribution function / direction-
                               *  specific densities of previous step.
                               */
    vava_val_t   v;        /*!<
                               * \brief Dimension-specific velocity (\f$
                               *  lu\,ts^{-1} \f$).
                               */
    va_val_t       dens;     //!< Density (\f$ mu\,lu^{-3} \f$).

    ///@}>

public:
    /*!
     * \brief Constructor.
     * \param prm - parameter file.
     */
    SCSP(const FilePrm &prm);

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
     *      loadRst;
     *      loadRst [ label="Load restart"];
     *    }
     *
     *    subgraph cluster1{
     *      label="Initialize";
     *      init -> svOpt1;
     *      init    [ shape=box, label="Apply initial values"];
     *      svOpt1  [ label="Save output file"];
     *    }
     *
     *    { svOpt1 loadRst } -> "t++";
     *
     *    subgraph cluster2{
     *      label="Main loop";
     *      subgraph cluster3{
     *        label="Core";
     *        style=filled;
     *        color=lightgray;
     *        node [ style=filled; color=white ];
     *        Streaming -> "Calculate macroscopic variables" -> BC -> Collision;
     *        BC [ label= "Apply boundary Conditions"];
     *      }
     *      "t++" -> Streaming;
     *      Collision -> isSaveOpt;
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
     * \brief Apply initial values (direction-specific densities).
     *
     * Initial direction-specfic densities are applied to void, gravity void,
     *  pressure boundary, velocity boundary, and bounceback boundary.
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
     * Applis to grids indexed "PB", "VB", "BB", and "SD".
     */
    void calcMacroVar();

    /*!
     * \brief The collision.
     *
     * Calculate \f$ f^{eq} \f$ in tems of \f$ u \f$, then perform the key step
     *  as described in the SCSP class description.
     *
     * The collosion step is applied to grid indexed "VD", "GV", "PB" and "VB".
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
     */
    void collision();

    /*!
     * \brief Apply Zou & He boundary conditions.
     *
     * The following boundary conditions are applied:
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
     *   {2\sum{\boldsymbol{f_{opp}}} + \sum{\boldsymbol{f_{normal}}}}
     *   {1-\boldsymbol{u}} \f]
     *  - Calculation of unknown distribution function:
     *   \f[ f_{i} = f_{opp} + 6w_i\rho\frac{\boldsymbol{u}}{c} +
     *   \frac{1}{2}\sum{f_{normal}e_i}\f]
     *   \f$ f_{opp} \f$ - opposit direction \n
     *   \f$ f_{normal} \f$ - directions perpendicular to the flow direction.
     */
    void zouHeBoundary();
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
     * \param dens - density.
     * \param v - velocity.
     * \return true if satisfied ending condition.
     */
    bool report(std::ostream &os,
                  const va_val_t& dens, const vava_val_t &v);

    ///@}

    /// \name Equilibrium functions.
    ///@{
    /*!
     * \brief Apply single rexlation time (SRT) equilibrium, i.e.
     *  Bhatnagar-Gross-Krook (BGK) approximation.
     *
     * The SRT model implements the following streaming and collision scheme:
     *  \f[\underbrace{f_i(\boldsymbol{x} + \boldsymbol{e}_i\delta t, t +
     *   \delta t) = f_i(\boldsymbol{x},t)}_{Streaming}
     *      - \underbrace{\frac{f_i(\boldsymbol{x},t)
     *                     - f_i^{eq}(\boldsymbol{x},t)}{\tau}}_{Collision} \f]
     *  The equilibrium distribution function \f$ f^{eq} \f$ is defined as: \f[
     *  f_i^{eq}(\boldsymbol{x}) = w_i\rho(\boldsymbol{x})
     *  \left[1+3\frac{\boldsymbol{e}_i\cdot\boldsymbol{u}}{c^2}
     *        + \frac{9}{2}\frac{(\boldsymbol{e}_i\cdot\boldsymbol{u})^2}{c^4}
     *        - \frac{3}{2}\frac{\boldsymbol{u}^2}{c^2}
     *  \right] \f]
     */
    void equilSRT();

    /*!
     * \brief Apply multiple rexlation time equilibrium (MRT) for D2Q9 model.
     *
     * The equilibrium function for MRT is:
     *  \f[\boldsymbol{f}(\boldsymbol{x}+\boldsymbol{e}_i\delta t, t + \delta t)
     *    - \boldsymbol{f}(\boldsymbol{x},t)
     *      = S\left[\boldsymbol{f_{eq}}(\boldsymbol{x},t)
     *             - \boldsymbol{f}(\boldsymbol{x},t)\right]
     *      = M^{-1}\hat{S}\left[M\boldsymbol{f^{eq}}(\boldsymbol{x},t)
     *        - M\boldsymbol{f}(\boldsymbol{x},t)\right] \f]
     *  where: \f$ S \f$ is a collision matrix, \f$ M \f$ is a transformation
     *  matrix, \f$ \hat{S} = M \cdot S \cdot M^{-1}\f$. \n
     * The row vectors in transformation matrix
     *  \f$ M \f$ is orthogonal, i.e. \f$ M\cdot M^{T} \f$ is a diagonal matrix.
     *  So the inverse of \f$ M^{-1} = M^T \cdot (M \cdot M^T)^{-1} \f$.
     *
     * For D2Q9:
     *  \f[ \boldsymbol{\hat{S}} = diag(1, s_{\nu}, s_{\nu}, 1, s_q, 1, s_q,
     *   s_{\nu}, s_{\nu}) \f]
     *  \f[ \boldsymbol{m^{eq}} = M\boldsymbol{f^{eq}} = \begin {bmatrix}
     *   \rho & -2\rho+3\rho_0(u_x^2 + u_y^2) & \rho - 3\rho_0(u_x^2+u_y^2) &
     *   \rho_0 u_x & -\rho_0 u_x & \rho_0 u_y & -\rho_0 u_y &
     *   \rho_0(u_x^2 - u_y^2) & \rho_0 u_x u_y
     *  \end {bmatrix} \f]
     *
     * For D3Q19:
     * \f[ \boldsymbol{\hat{S}}
     *  = diag(0, s_{e}, s_{\varepsilon}, 0, s_q, 0, s_q,
     *  0, s_q, s_{\nu}, s_{\pi}, s_{\nu}, s_{\pi}, s_{\nu}, s_{\nu}, s_{\nu},
     *  s_m, s_m, s_m) \\
     *  = diag(0, s_{\nu}, s_{\nu}, 0, s_q, 0, s_q, 0, s_q,
     *  s_{\nu}, s_{\nu}, s_{\nu}, s_{\nu}, s_{\nu}, s_{\nu}, s_{\nu},
     *  s_q, s_q, s_q)  \f]
     *  where \f$s_e = s_{\varepsilon} = s_{\pi} = s_{\nu} = \frac{1}{\tau},
     *  s_m = s_q]= 8 \frac{2-s_{\nu}}{8-s_{\nu}} \f$.
     *  \f[ \boldsymbol{m^{eq}} = M\boldsymbol{f^{eq}} = \begin {bmatrix}
     *  \rho & -11\rho+19\rho_0(u_x^2 + u_y^2 + u_z^2) & \frac{1}{2}(6\rho -
     *  11\rho_0(u_x^2 + u_y^2 + u_z^2)) & \rho_0 u_x & - \frac{2\rho_0 u_x}{3}
     *  & \rho_0 u_y & - \frac{2\rho_0 u_y}{3} & \rho_0 u_z &
     *  - \frac{2\rho_0 u_z}{3} & \rho_0(2u_x^2 - u_y^2 - u_z^2) &
     *  \frac{1}{2}\rho_0(-2u_x^2 + u_y^2 + u_z^2) & \rho_0(u_y^2 - u_x^2) &
     *  \frac{1}{2}\rho_0(- u_y^2 + u_z^2) & \rho_0 u_x u_y & \rho_0 u_y u_z &
     *  \rho_0 u_x u_z & 0 & 0 & 0
     *  \end {bmatrix} \f]
     */
    void equilMRT();

    ///@}
};

} // namespace LBM

} // namespace DRL

#endif // SCSP_H
