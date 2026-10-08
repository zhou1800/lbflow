/*! \file Bera-122x102y102z.h
 *
 * A parameter file of Bera sandstone with dimension of 102z102y122x. This file
 *  is for testing, and to demonstrat the parameter file.
 *
 * \warning This file is not a c/c++ header file. The extension (.h) is
 *  intended to highlight in vim and to parse by doxygen.
 *
 * \author Yu Zhou
 */

/*!
 * \page para Parameter description file.
 *
 * \brief The parameter description for input file.
 *
 * \note
 * 1. Keywords should be unique, otherwise, the later keyword overrides the
 *    previous one.
 * 2. Keywords should start with an alphabetic or a number, otherwise it is
 *    ignored.
 * 3. Keywords of pure integer is considered as index.
 * 4. The line following keyword is considered as values. The value type should
 *    be consistent in a line. Possible value types are: integer, string
 *    (enclosed with ""), and float.
 * 5. Multiple values (array) should be keep in a single line and separated by
 *    space.
 *
 * \bug When parsing string value, everything enclosed with "" will be loaded,
 *  for example "something" after ; and // will also be load.
 *
 * \note Parameters not required will be ignored by program.
 */

/*!
 * # Basic #
 * ## Grid ##
 * The grid defination.
 *
 * \param gridFile - grid file name (_without path, and working directory is
 *  "./" by default._).
 * \param dimension - grid dimension, in Z, Y, X order.
 * \param directions - directions, i.e. 19 for D3Q19, and 9 for D2Q9.
 * \param fileType - grid file type, i.e. "uint32", "uint16, and "uint8".
 * \param latticeLen - distance between two lattices (unite m).
 *
 ~~~~~~~~~~~~~~
    #grid
    {
        gridFile     "Bera-102z102y122x.raw";
        dimension    102 102 122;
        directions   19;
        fileType     "uint8";
        latticeLen   0.000001;
    }
 ~~~~~~~~~~~~~~
 *
 * ## Simulation control ##
 * These parameters controls the simulator.
 *
 * \param startTime - (ts) start time step for the simulation. (_While
 *  restarting, startTime should be equal to the number of old file, and the old
 *  file will not be overwritten._)
 * \param endTime - (ts) end time step for the simulation.
 * \param writeInterval - (ts) write output file on time step divisle by
 *  writeInterval.
 * \param reportInterval - (ts) report on time step divisle by reportInterval.
 * \param restartInterval - (ts) write restart file on time step divisle by
 *  restartInterval.
 * \param stopError - threshold error for stopping the simulation.
 * \param outPath - path and file name prefix of output file.
 * \param logPath - path and file name prefix of the log file.
 * \param reportDir - normal direction for report slice.
 * \param reportSlice - position for report slice.
 *
 ~~~~~~~~~~~~~~
    #control
    {
        startTime       0;
        endTime         200;
        writeInterval   10;
        reportInterval  10;
        restartInterval 10;
        stopError       0.000001;
        outPath         "results/Bera-122x102y102z";
        logPath         ".";
        reportDir       1;
        reportSlice     12 110;
    }
 ~~~~~~~~~~~~~~
 *
 * # Fluid model #
 * Avaliable models:
 * - scsp: single component, single phase.
 * - color: color fluid, two phases, two component
 * - shanchen: multiple component, multiple phase.
 *
 * ## Single component single phase model ##
 *
 * \param fluidModel - name of the fluid model.
   \param equilModel - choice of the equilibrium function. SRT or MRT.
 * \param rexlationTime - (\f$ \tau (ts) \f$) Relaxation time. Calculated by
 *  Kinematic Viscosity \f$ \nu (lu^2\,ts^{-1}) \f$. Given by: \f$
 *  \nu=c_s^2(\tau-\frac{1}{2}) \Rightarrow \tau = \frac{\nu}{c_s^2} +
 *  \frac{1}{2} \f$.
 * \param latticeSpeed - (\f$ c (lu/ts) \f$) lattice speed. Lattice speed \f$
 *  c = 1 \f$ and \f$ c_s = \sqrt{\frac{c}{3}}\,lu\,ts^{-1} \f$.
 * \param density - (\f$ \rho (mu\,lu^{-3}) \f$) density of fluid.
 *
 ~~~~~~~~~~~~~~
    #scsp
    {
        //fluidModel          "scsp";
        //equilModel          "SRT";
        //rexlationTime       1;
        //latticeSpeed        1;
        //density             1;
    }
 ~~~~~~~~~~~~~~
 *
 * ## Color fluid model property ##
 *
 * \param fluidModel - name of the fluid model.
   \param equilModel - choice of the equilibrium function. SRT or MRT.
 * \param rexlationTime - (\f$ \tau (ts) \f$) Relaxation time. Calculated by
 *  Kinematic Viscosity \f$ \nu^k (lu^2\,ts^{-1}) \f$. Given by: \f$
 *  \nu^k=\frac{c^2}{3}(\tau^k - \frac{1}{2})\delta_t \Rightarrow \tau^k =
 *  \frac{3 \nu^k}{c^2 \tau^k} + \frac{1}{2} \f$.
 * \param latticeSpeed - (\f$ c (lu/ts) \f$) lattice speed. \f$ c = 1 \f$
 *  For D2Q9 model: \f$ (c_s^k)^2 = \frac{3(1-\alpha_k)}{5}c \f$;
 *  for D3Q19 model: \f$ (c_s^k)^2 = \frac{1 - \alpha_k}{2}c \f$.
 * \param alpha0 - free parameter determins the sound speed for phase 0.
 *  \f$ \lambda = \frac{\rho_R}{\rho_B} = \frac{1-\alpha_B}{1- \alpha_R} \f$.
 * \param density - (\f$ mu/lu^3 \f$) density for each phase.
 * \param contactAngle - (\f$ \theta \f$) contact angle of phase 1.
 * \param surfaceTension - surface tension.
 * \param interfaceThickness - interface thickness.
 *
 ~~~~~~~~~~~~~~
    #colorfluid
    {
        fluidModel          "color";
        equilModel          "MRT";
        rexlationTime       0.57    1.215416667;
        latticeSpeed        1;
        alpha0              0.33333333333333333333;
        density             1       0.831683168;
        contactAngle        150;
        surfaceTension      0.085920139;
        interfaceThickness  0.7;
    }
 ~~~~~~~~~~~~~~
 *
 * ## S/C model property ##
 *
 * \param fluidModel - name of the fluid model.
 * \param components - name of components, i.e. "water", "oil", "gas".
 *
 ~~~~~~~~~~~~~~
    #SC
    {
        //fluidModel          "shanchen"
        //components  "oil"   "water";
    }
 ~~~~~~~~~~~~~~
 *
 * # Grid index #
 * Defination of indexes in grid file.
 *
 * ## Index values ##
 * - index: only numbers accepted, the number corresponding with values in grid
 *   file.
 * - Boundary type (type):
 *      - "VD" for void;
 *      - "GV" for void with gravity incorporated;
 *      - "PB" for pressure boundary;
 *      - "VB" for velocity boundary;
 *      - "AB" for axisymmetric boundary; (unfinished)
 *      - "BB" for bounce back boundary; (deprecated)
 *      - "HF" for half permeable film, 'Phase' for permeable phase;
 *      - "SD" for solid.
 * - Value of force (value):
 *      - force value (float).
 * - Force direction (dir):
 *      - Index of force direction (DmQn scheme),
 *        i.e. 1:x+; 2:x-; 3:y+; 4:y-; 5:z+; 6:z-.
 *      - For boundary, force direction must point to the inside for deducting
 *        unknown variables. e.g. at east boundary, dir should pointing to west.
 *      - For gravity, direction could be dispersed to every direction.
 * - Phase (phase):
 *      - Index of phase.
 *
 * ## Values for each type ##
 * | Type | value | dir | phase |
 * |------|-------|-----|-------|
 * |  VD  |       |     |   Y   |
 * |  GV  |   Y   |  Y  |   Y   |
 * |  PB  |   Y   |  Y  |   Y   |
 * |  VB  |   Y   |  Y  |   Y   |
 * |  AB  |       |  Y  |       | (unfinished)
 * |  BB  |       |     |       | (deprecated)
 * |  HF  |       |     |   Y   |
 * |  SD  |       |     |       |
 *
 * \note Periodic boundary condition is automatically applied. Specify a
 *  boundary to avoid periodicity.
 *
 ~~~~~~~~~~~~~~
    #Index
    {
      //        0       1             2       3
      //index   type    value         dir     phase
        0       "VD"    "N/A"         "N/A"   "0"
        50      "VB"    "0.000589167" "1"     "1"
        90      "VD"    "N/A"         "N/A"   "1"
        150     "PB"    "1"           "2"     "0"
        128     "BB"    "N/A"         "N/A"   "N/A"
        255     "SD"    "N/A"         "N/A"   "N/A"
    }
 ~~~~~~~~~~~~~~
 */
