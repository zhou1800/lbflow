/*! \file fileprm.h
 * \author Yu Zhou
 */
#ifndef FILEPRM_H
#define FILEPRM_H

#include "file.h"
#include <stdexcept>

namespace DRL {

namespace IO {

/*!
 * \brief The FilePrm class is a parallel parameter file parser.
 *
 * The parameter file is constructed with lines of key-value pair. Key should be
 *  unique, started with alpha or number. Value could be a string, integer,
 *  float, and array could either be array of number or string.
 *
 * Fratures:
 *  - Load parameters by master core, then sync between cores.
 *  - Check parameter type and name while compiling
 *  - Runtime error thrown on required parameter missing in file. (by using
 *    parameter pointers)
 *
 * \note For adding new parameter, follow these steps:
 *  1. Add a new private variable and coorsponding get function.
 *  2. Update the new variable in clear(), parseAll(), and Constructors'
 *     initialization list.
 *  3. Done.
 */
class FilePrm : public File
{
public:
    /// \name Index definations.
    ///@{
    /*!
     * \brief The IDXTYPE enum defines the acceptable index type.
     *
     * The enum value of index types are assigned by \f$ 2^n \f$ to get a unique
     *  sum of any type combinations.
     */
    enum IDXTYPE {
        VD  = 0x1,  //!< Void, "VD" in parameter file.
        GV  = 0x2,  //!< Void with gravity incorported, "GV" in parameter file.
        PB  = 0x4,  //!< Pressure boundary, "PB" in parameter file.
        VB  = 0x8,  //!< Velocity boundary, "VB" in parameter file.
        AB  = 0x10, //!< Axisymmatric bouhndary, "AB" in parameter file.
        BB  = 0x20, //!< Bounce back boundary, "BB" in parameter file.
        HF  = 0x40, //!< Half-permeable film, "HF" in parameter file.
        SD  = 0x80, //!< Solid, "SD" in parameter file.
    };

    /*!
     * \brief The IDXVAL struct defines the content value of a index.
     *
     * A index is a value of pixles in the geomotery file. Index value is
     *  attributes related to that index. These attributes are defined in
     *  parameter file, i.e. 0 represents void, and 255 represents matrix.
     */
    struct IDXVAL{
        IDXTYPE Type;           //!< Index type.
        double  Value;          //!< Value of force.
        size_t  Dir;            //!< Force direction.
        size_t  Phase;          //!< Phase.
        /*!
         * \brief Constructor.
         * \param val - string array of value in the same order as in the struct
         *  defination.
         */
        IDXVAL(const vector<string> & val);
    private:
        static constexpr size_t N_IDXPARM {4};  //!< Number of parameters.
    };
    ///@}

private:
    map<string, string> iPrm;           //!< Parameter index, pair of key-value.

    /// \name Parameters in input file: #grid mesh
    ///@{
    string*         gridFile;       //!< Grid filename.
    va_idx_t        dimension;      //!< Grid dimension.
    size_t*         directions;     //!< Number of directions.
    string*         fileType;       //!< Grid file type.
    double*         latticeLen;     //!< Lattice length (m).

    ///@}

    /// \name Parameters in input file: #simulation control
    ///@{
    size_t*         tStart;         //!< Start time.
    size_t*         tEnd;           //!< End time.
    size_t*         tStep;          //!< Time step for output result file.
    size_t*         tRpt;           //!< Time step for report.
    size_t*         tRst;           //!< Time step for output restart file.
    double*         stopError;      //!< Stop criterion of minium error.
    string*         outPath;        //!< Output file path.
    string*         logPath;        //!< Log file path.
    size_t*         reportDir;      //!< Normal direction for report slice.
    vector<size_t>  reportSlice;    //!< Position for report slice.

    ///@}

    /// \name Parameters in input file: #property (scsp)
    ///@{
    string*         fluidModel;     //!< Fluid model.
    string*         equilModel;     //!< Equilibrum Model, i.e. MRT or SRT.
    vector<double>  rexlationTime;  //!< Rexlation time.
    double*         latticeSpeed;   //!< Lattice speed \f$c(lu\,ts^{-1})\f$.
    double*         alpha0;         //!< \f$ \alpha[0] \f$
    vector<double>  density;        //!< Density.

    ///@}

    /// \name Parameters in input file: #property (color-fluid)
    ///@{
    vector<string>   phase;         //!< Name of phases.
    double*          cAngle;        //!< Contact angle.
    double*          sTension;      //!< Surface tension.
    double*          iThickness;    //!< Interface thickness.

    ///@}

    /// \name Parameters in input file: #property (S/C)
    ///@{
    vector<string>   components;    //!< Name of components.

    ///@}

    /// \name Parameters in input file: #index
    ///@{
    using idx_set = set<size_t>;        //!< Set of indexes.
    map<size_t, IDXVAL*>  idx;          //!< Indexes and contents.
    mutable
    map<IDXTYPE, idx_set> idxByType;     /*!<
                                         * \brief Index lookup table by type.
                                         *
                                         * For fast access to indexes with the
                                         *  same index type.
                                         */
    mutable
    map<size_t, idx_set>  idxByTypComb; /*!<
                                         * \brief Index lookup table by type
                                         *  combination.
                                         *
                                         * This attribute is updated on the fly.
                                         *
                                         * \sa getIdxByType() to get and update
                                         *  this parameter.
                                         */

    ///@}

public:
    /// \name Constructor and destructor
    ///@{
    /*!
     * \brief Construct with a parameter file name for loading.
     *
     * All parameter pointers are initialized to 0.
     *
     * \note
     * By setting a parameter pointer to 0 indicates that this parameter
     *  has not been or failed to be load in parameter file.
     *
     * \note To add a new parameter, add it to this constructor's initialization
     *  list. (step 1/4)
     *
     * \param filename - parameter file name
     */
    FilePrm(const string& filename);

    /*!
     * \brief Destructor.
     *
     * Delete all parameter pointers.
     *
     * \sa Clean()
     */
    virtual ~FilePrm();
    ///@}

    /// \name Parameter file manipulation
    ///@{
    /*!
     * \brief Open a file for loading parameters.
     *
     * More than a IO::Open(), it also reset all parameters to 0.
     *
     * \sa clear() the old parameters.
     * \param filename - parameter file name.
     */
    virtual void Open(const string& filename) override;

    // Validate simulation parameters before constructing or indexing a grid.
    void ValidateSimulation() const;
    string getWallGradient() const;
    string getViscosityInterpolation() const;
    string getCapillaryModel() const;
    string getWallRecoloring() const;
    string getCurvatureStencil() const;
    string getPhaseBoundary() const;
    string getRecoloringModel() const;
    string getMrtSpectrum() const;

    /// \name Get functions for #grid mesh
    ///@{
    /*!
     * \brief Get parameter: grid file name (gridFile).
     * \return gride file name.
     */
    const string& getGFile() const { return get(gridFile); }

    /*!
     * \brief Get parameter: length of all dimension (dimension).
     * \return length of all dimension.
     */
    const va_idx_t& getD() const { return dimension; }

    /*!
     * \brief Get parameter: number of directions (directions).
     * \return number of directions.
     */
    size_t getQ() const { return get(directions); }

    /*!
     * \brief Get file type.
     * \return grid file type.
     */
    const string& getFileType() const { return get(fileType); }

    /*!
     * \brief Get lattice length.
     * \return lattice length (m).
     */
    double getLatticeLen() const { return get(latticeLen); }
    ///@}

    /// \name Get functions for #simulation control
    ///@{
    /*!
     * \brief Get parameter: start time step (tStart).
     * \return start time step.
     */
    size_t getTStart() const { return get(tStart); }

    /*!
     * \brief Get parameter: end time step (tEnd).
     * \return end time step
     */
    size_t getTEnd() const { return get(tEnd); }

    /*!
     * \brief Get parameter: result file output time step (tStep).
     * \return differ time step.
     */
    size_t getTStep() const { return get(tStep); }

    /*!
     * \brief Get parameter: report time step (tStep).
     * \return differ time step.
     */
    size_t getTRpt() const { return get(tRpt); }

    /*!
     * \brief Get parameter: restart file output time step (tRst).
     * \return differ time step.
     */
    size_t getTRst() const { return get(tRst); }

    /*!
     * \brief Get stop criterion of minuum error.
     * \return stop error.
     */
    double getStopError() const { return get(stopError); }

    /*!
     * \brief Get parameter: output path (outPath).
     * \return output path.
     */
    const string& getOutPath() const { return get(outPath); }

    /*!
     * \brief Get log file path (and file name).
     * \return log file path.
     */
    const string& getLogPath() const { return get(logPath); }

    /*!
     * \brief Get normal direction for report slice.
     * \return normal direction for report slice.
     */
    size_t getReportDir() const { return get(reportDir); }

    /*!
     * \brief Get position for report slice.
     * \return position for report slice.
     */
    const vector<size_t>& getReportSlice() const { return reportSlice; }


    ///@}

    /// \name Get functions for #property (scsp)
    ///@{
    /*!
     * \brief Get parameter: name of fluid model (fluidModel).
     * \return name of fluid model.
     */
    const string& getFluidModel() const { return get(fluidModel); }

    const string& getEquilModel() const { return get(equilModel); }

    /*!
     * \brief Get parameter: rexlation time (rexlationTime).
     * \return relxation time of ith phase.
     */
    const vector<double>& getRTime() const { return rexlationTime; }

    /*!
     * \brief Get parameter: lattice speed (latticeSpeed).
     * \return lattice speed.
     */
    double getLSpeed() const { return get(latticeSpeed); }

    /*!
     * \brief Get \f$ \alpha_0 \f$ for colorfluid model, which determins the
     *  sound of speed.
     *
     * \return the free parameter \f$ \alpha_0 \f$ which determins the sound
     *  speed.
     */
    double getAlpha0() const { return get(alpha0); }

    /*!
     * \brief Get parameter: density (density).
     * \return density of ith phase.
     */
    const vector<double>& getDensity() const { return density; }
    ///@}

    /// \name Get functions for #property (color-fluid)
    ///@{
    /*!
     * \brief Get parameter: name of all phases (phase).
     * \return names of all phases.
     */
    const vector<string>& getPhase() const { return phase; }

    /*!
     * \brief Get parameter: contact angle (cAngle).
     * \return contact angle.
     */
    double getCAngle() const { return get(cAngle); }

    /*!
     * \brief Get parameter: surface tension (sTension).
     * \return surface tension.
     */
    double getSTension() const { return get(sTension); }

    /*!
     * \brief Get parameter: interface thickness (iThickness).
     * \return interface thickness.
     */
    double getIThickness() const { return get(iThickness); }
    ///@}

    /// \name Get functions for #property (S/C)
    ///@{
    /*!
     * \brief Get parameter: name of all components (components).
     * \return names of all components
     */
    const vector<string>& getComponents() const { return components; }
    ///@}

    /// \name Get functions for #index
    ///@{
    /*!
     * \brief Get index content.
     *
     * \return index content.
     */
    const IDXVAL& getIVal(size_t i) const { return *(idx.at(i)); }

    /*!
     * \brief Get the color indexes array by index types.
     *
     * Return an empty array if all given types has no indexes.
     *
     * First try to get indexes from pre-caculated attibute idxByTypComb, if
     *  failed, build the indexes set for given type combination by idxByType,
     *  and store it in idxByTypComb for later use.
     *
     * \param types - unique index type combinations.
     * \return indexes of target types.
     */
    const set<size_t> &getIdx(const set<IDXTYPE>& types) const;
    ///@}

private:
    /// \name Helper functions for loading a new parameter file.
    ///@{
    /*!
     * \brief Load parameters.
     *
     * Following functions has been performed to load the parameters:
     *  1. clear():     Reset all paremeter pointers and parameter index.
     *  2. loadFile():  Load file into parameter index (iPrm) by core 0
     *  3. sync():      Distribute parameter index to each core.
     *  4. parseAll():  Parse parameters on all cores. Store within the class.
     *  5. Done!
     *
     * \sa clear(), loadFile(), sync(), parseAll().
     */
    void load();

    /*!
     * \brief Load parameter file content into parameter index(iPrm).
     *
     * Read key-value pair. Every line start with alpha or number is considered
     *  as a key. Keys are tied with value, which is a following string extend
     *  to the end of this line.
     *
     * \exception runtime error on failed to load the file.
     */
    void loadFile();

    /*!
     * \brief Syncronize parameter index iPrm between cores.
     *
     * Syncronizs following these steps:
     * 1. Transfer iPrm into strings for syncronize.
     * 2. Broad cast number of strings, length of each string.
     * 3. Prapare mpi packed structure of parameter strings for sending and
     *    receiving using MPI.
     * 4. Syncronizs strings.
     * 5. Transfer strings back to iPrm.
     * 6. Clean up.
     */
    void sync();

    /*!
     * \brief Parse the parameter index into parameters.
     *
     * The intention of designing a pre-parsed mechanism before useage (instead
     *  of a runtime parser which accept string of parameter name as a input
     *  parameter and search for the value), is to provent runtime error from
     *  happing (caused by careless coding).
     *
     * When accessing an parameter, the following mistakes will be noticed
     *  during compling.
     *
     * - wrong parameter name,
     * - wrong parameter type,
     * - undifined parameter.
     *
     * \note To add a new parameter, add the Key to this function manually.
     *  (step 2/4)
     *
     * \sa parse(const string& str, T **val)const
     */
    void parseAll();

    /*!
     * \brief Clean up all parameters.
     *
     * Delete all parameter pointers, and set all the pointers to 0.
     *
     * \note
     * By setting a parameter pointer to 0 indicates that this parameter
     *  has not been load or failed to be load.
     *
     * \note To add a new parameter, add the delete opertion in this function
     *  manually. (step 3/4)
     */
    void clear();
    ///@}

    /// \name Helper functions for getting value.
    ///@{
    template <typename T>
    const T& get(T* val) const;

    ///@}

    /// \name Parameters' parser (helper functions for parseAll())
    ///@{
    /*!
     * \brief Parse Entrance.
     *
     * This function allocate memory for value, then call parse().
     * Automatically delete the value and set pointer to 0
     * on failed parse.
     *
     * \sa parse(const string& str, size_t* i)const
     * \sa parse(const string& str, double* d)const
     * \sa parse(const string& str, string* s)const
     * \sa parse(const string& str, vector<size_t>* ai)const
     * \sa parse(const string& str, vector<double>* ad)const
     * \sa parse(const string& str, vector<string>* as)const
     * \param str - string to be parsed.
     * \param val - value used to store the parsed result.
     */
    template <typename T> static
    void parse(const string& str, T **val);

    /*!
     * \brief Parse an integer.
     * \param str - string to be convert.
     * \param[out] i - result integer.
     * \return true if parse is successful.
     */
    static bool parse(const string& str, size_t *i);

    /*!
     * \brief Parse an double percision float.
     * \param str - string to be convert.
     * \param[out] d - result float.
     * \return true if parse is successful.
     */
    static bool parse(const string& str, double* d);

    /*!
     * \brief Parse a string.
     *
     * String is enclosed with "".
     *
     * \param str - string to be convert.
     * \param[out] s result string.
     * \return true if parse is successful.
     */
    static bool parse(const string& str, string* s);

    /*!
     * \brief Parse an array of integer.
     * \param str - string to be convert.
     * \param[out] ai - result array of integer.
     * \return true if parse is successful.
     */
    static bool parse(const string& str, vector<size_t>* ai);

    /*!
     * \brief Parse an array of double precision float.
     * \param str - string to be convert.
     * \param[out] ad - result array of float.
     * \return true if parse is successful.
     */
    static bool parse(const string& str, vector<double>* ad);

    /*!
     * \brief Parse an array of string.
     *
     * \bug Anything enclosed within "" will be load, * i.e.
     *  "something" after ; and // will also be load.
     *
     * \param str - string to be convert.
     * \param[out] as - result array of string.
     * \return true if parse is successful.
     */
    static bool parse(const string& str, vector<string>* as);

    /*!
     * \brief Try to parse a pair of string into index.
     *
     * Convert the pair of string to index headers and add them to the index
     *  head map. It also update the index lookup table by type (idxByType),
     *  which is a  Look up table for fast access to all color indexes of the
     *  same index type.
     *
     * \param a - a line of parameter.
     */
    void parse(const std::pair<string, string>& a);
    ///@}

    /// \name Helper function for Sync()
    ///@{
    // transfer between PrmIndex(map<string, string>) and chars
    /*!
     * \brief Convert parameter index to strings.
     *
     * Convert to a sequence of key - value strings. Also log the length
     * of string and number of strings. This function is inverse of Str2Prm.
     *
     * \sa str2prm()
     * \param[out] strs - content of strings.
     * \param[out] strLen - lengths of strings
     * \param[out] strCnt - number of string
     */
    void prm2str(char ***strs, int **strLen, unsigned* strCnt) const;

    /*!
     * \brief Convert strings to parameter index.
     *
     * Convert strings to parameter index. Strings must be end with \0.
     * This function is inverse of Prm2Str().
     *
     * \sa prm2str()
     * \param strs - Strings to be convert.
     * \param strCnt - number of strings.
     */
    void str2prm(char** strs, unsigned strCnt);
    ///@}

    /// \name Helper functions for parse index.
    ///@{
    /*!
     * \brief Convert string to double percision float value.
     *
     * A specified version for getting double value from parameter string.
     *  - With convert successful confirmation.
     *  - With dafault value returning.
     *
     * \exception invalid_argument on failed to parse string and isMust is true.
     * \param val - string to be converted.
     * \param isMust -whether to throw exception on invalid string of double.
     * \return double presision value (_return 0(as default value) if string is
     *  invalid and iMust == false_).
     */
    static
    double str2d(const string& val, bool isMust = true);

    /*!
     * \brief Convert string to unsigned long.
     *
     * A specified version for getting unsigned long from parameter string.
     *  - With convert successful confirmation.
     *  - With dafault value returning.
     *
     * \exception invalid_argument on failed to parse string and isMust is true.
     * \param val - string to be converted.
     * \param isMust - whether to throw exception on invalid string of unsigned
     *  long.
     * \return unsigned long value (_return 0(as default value) if string is
     *  invalid and iMust == false_).
     */
    static
    size_t str2ul(const string& val, bool isMust = true);
    ///@}
};

using IT = FilePrm::IDXTYPE;   //!< Grid index type.

template <typename T> inline
const T &FilePrm::get(T *val) const
{
    if (!val)
        throw std::runtime_error("the value may not be defined in parameter file.");
    return *val;
}

template <typename T> inline
void FilePrm::parse(const string &str, T **val)
{
    if (!*val)
        *val = new T;
    if (!parse(str, *val)) {
        delete *val;
        *val = 0;
    }
}

} // namespace IO

} // namespace DRL

#endif // FILEPRM_H
