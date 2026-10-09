/*! \file fileprm.cpp
 * \author Yu Zhou
 */
#include "fileprm.h"
#include "thread.h"
#include <fstream>
#include <sstream>
#include <cstring>
#include <numeric>
#include <cmath>
#include <limits>
#include <algorithm>

namespace DRL {

namespace IO {

using namespace std;

FilePrm::FilePrm(const string &filename)
    : File(filename), iPrm()
    , gridFile(0), dimension(), directions(0), fileType(0), latticeLen(0)
    , tStart(0), tEnd(0), tStep(0), tRpt(0), tRst(0), stopError(0)
    , outPath(0), logPath(0), reportDir(0), reportSlice()
    , fluidModel(0), equilModel(0)
    , rexlationTime(), latticeSpeed(0), alpha0(0), density()
    , phase(), cAngle(0), sTension(0), iThickness(0)
    , components()
    , idx()
{
    load();
}

FilePrm::~FilePrm()
{
    clear();
}

void FilePrm::Open(const string &filename)
{
    File::Open(filename);
    clear();
    load();
}

string FilePrm::getWallGradient() const
{
    const auto value = iPrm.find("wallGradient");
    if (value == iPrm.end()) return "legacy";
    string result;
    if (!parse(value->second, &result))
        throw invalid_argument("wallGradient must be a quoted string.");
    return result;
}

string FilePrm::getViscosityInterpolation() const
{
    const auto value = iPrm.find("viscosityInterpolation");
    if (value == iPrm.end()) return "arithmetic";
    string result;
    if (!parse(value->second, &result))
        throw invalid_argument("viscosityInterpolation must be a quoted string.");
    return result;
}

string FilePrm::getCapillaryModel() const
{
    const auto value = iPrm.find("capillaryModel");
    if (value == iPrm.end()) return "perturbation";
    string result;
    if (!parse(value->second, &result))
        throw invalid_argument("capillaryModel must be a quoted string.");
    return result;
}

void FilePrm::ValidateSimulation() const
{
    const size_t ndim = dimension.size();
    const size_t q = getQ();
    if (!((ndim == 2 && q == 9) || (ndim == 3 && q == 19)))
        throw invalid_argument("Only D2Q9 and D3Q19 are validated for simulation.");
    size_t count = 1;
    for (size_t n : dimension) {
        if (n == 0 || n > numeric_limits<size_t>::max() / count)
            throw invalid_argument("Invalid grid dimensions.");
        count *= n;
    }
    if (getFileType() != "uint8")
        throw invalid_argument("Simulation geometry must use uint8 labels.");
    if (getTEnd() < getTStart() || getTEnd() == numeric_limits<size_t>::max() ||
        getTStep() == 0 || getTRpt() == 0 || getTRst() == 0)
        throw invalid_argument("Invalid time range or zero output/report/restart interval.");
    if (!isfinite(getStopError()) || getStopError() <= 0.0 ||
        !isfinite(getLatticeLen()) || getLatticeLen() <= 0.0)
        throw invalid_argument("Invalid convergence tolerance or lattice length.");
    if (getLSpeed() != 1.0)
        throw invalid_argument("This research implementation requires latticeSpeed=1.");
    const string& model = getFluidModel();
    if (model != "scsp" && model != "color")
        throw invalid_argument("Only scsp and color fluid models are implemented.");
    if (getEquilModel() != "SRT" && getEquilModel() != "MRT")
        throw invalid_argument("Equilibrium model must be SRT or MRT.");
    const size_t phases = model == "color" ? 2 : 1;
    if (rexlationTime.size() != phases || density.size() != phases)
        throw invalid_argument("Density and relaxation-time counts must match the fluid model.");
    for (size_t k = 0; k < phases; ++k)
        if (!isfinite(rexlationTime[k]) || rexlationTime[k] <= 0.5 ||
            !isfinite(density[k]) || density[k] <= 0.0)
            throw invalid_argument("Densities must be positive; relaxation times must exceed 0.5.");
    if (model == "color") {
        const string capillary = getCapillaryModel();
        if (capillary != "perturbation" && capillary != "csf")
            throw invalid_argument("capillaryModel must be perturbation or csf.");
        if (capillary == "csf") {
            const double standardAlpha = ndim == 2 ? 4.0/9.0 : 1.0/3.0;
            if (density[0] != density[1] || abs(getAlpha0()-standardAlpha) > 1e-14)
                throw invalid_argument("Experimental CSF requires equal reference densities and standard pressure weights.");
            for (const auto& item : idx) {
                const IT type = item.second->Type;
                if (type != IT::VD && type != IT::GV && type != IT::SD && type != IT::BB)
                    throw invalid_argument("Experimental CSF supports only closed VD/GV flow with SD/BB walls.");
            }
        }
        if (getViscosityInterpolation() != "arithmetic" && getViscosityInterpolation() != "harmonic")
            throw invalid_argument("viscosityInterpolation must be arithmetic or harmonic.");
        if (getWallGradient() != "legacy" && getWallGradient() != "renormalized" &&
            getWallGradient() != "tanh")
            throw invalid_argument("wallGradient must be legacy, renormalized or tanh.");
        const double a0 = getAlpha0();
        const double a1 = 1.0 - density[0] / density[1] * (1.0 - a0);
        if (!isfinite(a0) || a0 < 0.0 || a0 >= 1.0 ||
            !isfinite(a1) || a1 < 0.0 || a1 >= 1.0 ||
            !isfinite(getSTension()) || getSTension() < 0.0 ||
            !isfinite(getIThickness()) || getIThickness() < 0.0 || getIThickness() > 1.0)
            throw invalid_argument("Invalid phase weights, surface tension, or recoloring beta.");
    }
    if (getReportDir() >= q)
        throw invalid_argument("Report direction is outside the lattice.");
    if (getReportDir() != 0) {
        if (reportSlice.size() < 2)
            throw invalid_argument("At least two report slices are required.");
        // Validate against the selected axis after GridSpace maps the direction.
        if (*max_element(reportSlice.begin(), reportSlice.end()) >=
            *max_element(begin(dimension), end(dimension)))
            throw invalid_argument("Report slice is outside the grid.");
    }
    for (const auto& item : idx) {
        if (item.first > 255)
            throw invalid_argument("Geometry labels must fit in uint8.");
        const IDXVAL& v = *item.second;
        if (!isfinite(v.Value))
            throw invalid_argument("Non-finite boundary value.");
        if ((v.Type == IT::VD || v.Type == IT::GV || v.Type == IT::PB ||
             v.Type == IT::VB || v.Type == IT::HF) && v.Phase >= phases)
            throw invalid_argument("Boundary phase is outside the configured fluid model.");
        if ((v.Type == IT::GV || v.Type == IT::PB || v.Type == IT::VB) &&
            (v.Dir == 0 || v.Dir >= q))
            throw invalid_argument("Boundary direction must be a non-rest lattice direction.");
        if (v.Type == IT::PB && v.Value <= 0.0)
            throw invalid_argument("Pressure-boundary density must be positive.");
        if ((v.Type == IT::PB || v.Type == IT::VB) && v.Dir > 2 * ndim)
            throw invalid_argument("Pressure and velocity boundaries require an axial direction.");
        if (v.Type == IT::SD && (v.Value < 0.0 || v.Value > 180.0))
            throw invalid_argument("Solid contact angle must be between 0 and 180 degrees.");
        if (v.Type == IT::AB)
            throw invalid_argument("Axisymmetric AB boundary handling is not implemented.");
    }
}

void FilePrm::load()
{
    if (!isValid())
        throw runtime_error("FilePrm::Load() - Input file does not exist.");

    if (mpi().getId() == 0)     // for head core: load an input file
        loadFile();
    sync();
    parseAll();
}

void FilePrm::loadFile()
{
    // Load parameter index from file
    ifstream ifs(fname);
    if (!ifs.is_open())
        throw runtime_error("FilePrm::Load() - Failed to load input file");
    for (string key; ifs >> key;){
        string line;
        getline(ifs, line);
        if (isalnum(*key.begin()))
            iPrm[move(key)] = move(line);
    }
    ifs.close();
}

void FilePrm::sync()
{
    unsigned    strCnt;         // number of strings
    int*        strLen = 0;     // string length
    char**      strs = 0;       // content of each string

    // get strings from processor 0
    if (mpi().getId() == 0)
        prm2str(&strs, &strLen, &strCnt);

    // sync number of strings
    mpi().Bcast(&strCnt, 1, MPI_UNSIGNED, 0);

    // sync length of strings
    if (mpi().getId() != 0){
        strLen = new int[strCnt];
        strs = new char*[strCnt];
    }
    mpi().Bcast(strLen, strCnt, MPI_UNSIGNED, 0);

    // sync contents of strings
    if (mpi().getId() != 0){
        for(unsigned i = 0; i < strCnt; ++i)
            strs[i] = new char[strLen[i]];
    }
    MPI_Datatype strPacked;
    void* pFirst = mpi().Pack(&strPacked, strs, strLen, strCnt);
    mpi().Bcast(pFirst, 1, strPacked, 0);
    mpi().TypeFree(&strPacked);

    // put strings into processors
    if (mpi().getId() != 0)
        str2prm(strs, strCnt);

    // clean up
    for(size_t i = 0; i < strCnt; ++i)
        delete [] strs[i];
    delete [] strs;
    delete [] strLen;
}

void FilePrm::parseAll()
{
    for (const auto &a : iPrm){
        // #grid mesh
        if (a.first == "gridFile")
            parse(a.second, &gridFile);
        else if (a.first == "dimension"){
            vector<size_t> dim;
            parse(a.second, &dim);
            dimension.resize(dim.size());
            for (size_t i = 0; i < dim.size(); ++i)
                dimension[i] = dim[i];
        }
        else if (a.first == "directions")
            parse(a.second, &directions);
        else if (a.first == "fileType")
            parse(a.second, &fileType);
        else if (a.first == "latticeLen")
            parse(a.second, &latticeLen);
        // #simulation control
        else if (a.first == "startTime")
            parse(a.second, &tStart);
        else if (a.first == "endTime")
            parse(a.second, &tEnd);
        else if (a.first == "writeInterval")
            parse(a.second, &tStep);
        else if (a.first == "reportInterval")
            parse(a.second, &tRpt);
        else if (a.first == "restartInterval")
            parse(a.second, &tRst);
        else if (a.first == "stopError")
            parse(a.second, &stopError);
        else if (a.first == "outPath")
            parse(a.second, &outPath);
        else if (a.first == "logPath")
            parse(a.second, &logPath);
        else if (a.first == "reportDir")
            parse(a.second, &reportDir);
        else if (a.first == "reportSlice")
            parse(a.second, &reportSlice);
        // #property (scsp)
        else if (a.first == "fluidModel")
            parse(a.second, &fluidModel);
        else if (a.first == "equilModel")
            parse(a.second, &equilModel);
        else if (a.first == "rexlationTime")
            parse(a.second, &rexlationTime);
        else if (a.first == "latticeSpeed")
            parse(a.second, &latticeSpeed);
        else if (a.first == "alpha0")
            parse(a.second, &alpha0);
        else if (a.first == "density")
            parse(a.second, &density);
        // #property (color-fluid)
        else if (a.first == "phase")
            parse(a.second, &phase);
        else if (a.first == "contactAngle")
            parse(a.second, &cAngle);
        else if (a.first == "surfaceTension")
            parse(a.second, &sTension);
        else if (a.first == "wallGradient") {
            // Optional research selector is read directly from iPrm.
        }
        else if (a.first == "viscosityInterpolation") {
            // Optional interpolation selector is read directly from iPrm.
        }
        else if (a.first == "capillaryModel") {
            // Optional experimental capillary force selector.
        }
        else if (a.first == "interfaceThickness")
            parse(a.second, &iThickness);
        // #property (S/C)
        else if (a.first == "components")
            parse(a.second, &components);
        // #value index with its defination.
        else parse(a);
    }
}

void FilePrm::clear()
{
    iPrm.clear();
    delete gridFile;        gridFile = 0;
    dimension.resize(0);
    delete directions;      directions = 0;
    delete fileType;        fileType = 0;
    delete latticeLen;      latticeLen = 0;
    delete tStart;          tStart = 0;
    delete tEnd;            tEnd = 0;
    delete tStep;           tStep = 0;
    delete tRpt;            tRpt = 0;
    delete tRst;            tRst = 0;
    delete stopError;       stopError = 0;
    delete outPath;         outPath = 0;
    delete logPath;         logPath = 0;
    delete reportDir;       reportDir = 0;
    reportSlice.clear();
    delete fluidModel;      fluidModel = 0;
    delete equilModel;      equilModel = 0;
    rexlationTime.clear();
    delete latticeSpeed;    latticeSpeed = 0;
    delete alpha0;          alpha0 = 0;
    density.clear();
    phase.clear();
    delete cAngle;          cAngle = 0;
    delete sTension;        sTension = 0;
    delete iThickness;      iThickness = 0;
    components.clear();
    for (auto a : idx) delete a.second;
    idx.clear();
}

bool FilePrm::parse(const string &str, size_t *i)
{
    istringstream is{str};
    is >> *i;
    return !is.fail();
}

bool FilePrm::parse(const string &str, double *d)
{
    istringstream is{str};
    is >> *d;
    return !is.fail();
}

bool FilePrm::parse(const string &str, string *s)
{
    s->clear();
    istringstream is{str};
    char c;
    while (is.get(c) && c != '"');
    while(is.get(c) && c != '"')
        s->push_back(move(c));
    return !is.fail();
}

bool FilePrm::parse(const string &str, vector<size_t> *ai)
{
    ai->clear();
    istringstream is{str};
    while (is) {
        int i;
        is >> i;
        if (is.fail())
            break;
        ai->push_back(move(i));
    }
    return ai->size() != 0;
}

bool FilePrm::parse(const string &str, vector<double> *ad)
{
    ad->clear();
    istringstream is{str};
    while (is) {
        double i;
        is >> i;
        if (is.fail())
            break;
        ad->push_back(move(i));
    }
    return ad->size() != 0;
}

bool FilePrm::parse(const string &str, vector<string> *as)
{
    as->clear();
    istringstream is{str};
    while (is) {
        string s;
        char c;
        while (is.get(c) && c != '"');
        s.clear();
        while(is.get(c) && c != '"')
            s.push_back(move(c));

        if (is.fail())
            break;
        as->push_back(move(s));
    }
    return as->size() != 0;
}

void FilePrm::parse(const pair<string, string> &a)
{
    char* pEnd;
    size_t iVal = strtoul(a.first.c_str(), &pEnd, 10);
    if (pEnd == a.first.c_str())        // not a number
        return;
    vector<string>* iCont = 0;
    parse(a.second, &iCont);
    if (iCont) {
        idx[iVal] = new IDXVAL(*iCont);
        idxByType[idx[iVal]->Type].insert(iVal);
    }
}

void FilePrm::prm2str(char ***strs, int **strLen, unsigned *strCnt) const
{
    // get number of strings
    *strCnt = iPrm.size() * 2;

    // get length and contents of strings
    *strLen = new int[*strCnt];    // string length
    *strs = new char*[*strCnt];         // content of each string

    size_t iStr = 0;        // index of string
    if (mpi().getId() == 0) {
        for (const auto &pair : iPrm){
            for (const auto &a : {pair.first, pair.second}){
                size_t sz = a.length() + 1;
                (*strs)[iStr] = new char[sz];
                strcpy((*strs)[iStr], a.c_str());
                (*strLen)[iStr++] = sz;
            }
        }
    }
}

void FilePrm::str2prm(char **strs, unsigned strCnt)
{
    for (size_t i = 0; i < strCnt; i+=2)
        iPrm[move(strs[i])] = move(strs[i + 1]);
}

double FilePrm::str2d(const string &val, bool isMust)
{
    char* pEnd = 0;
    double value = strtod(val.c_str(), &pEnd);
    if (pEnd == val.c_str()) {
        if (isMust)
            throw invalid_argument("IdxHead::str2d() - "
                                   "double value not defined");
        else
            value = 0;   // return default value.
    }
    return value;
}

size_t FilePrm::str2ul(const string &val, bool isMust)
{
    char* pEnd = 0;
    size_t value = strtoul(val.c_str(), &pEnd, 10);
    if (pEnd == val.c_str()) {
        if (isMust)
            throw invalid_argument("IdxHead::str2ul() -"
                                   " unsigned long value not defined");
        else
            value = 0;  // return default value.
    }
    return value;
}

FilePrm::IDXVAL::IDXVAL(const vector<string> &val)
    : Type(IT::VD)
    , Value(0.0)
    , Dir(0)
    , Phase(0)
{
    if (val.size() < N_IDXPARM)
        throw invalid_argument("FilePrm::IDX::IDX() - insufficient val size.");

    if (val[0] == "VD"){        Type = IT::VD;
        Phase = str2ul(val[3]);
    }
    else if (val[0] == "GV"){   Type = IT::GV;
        Value = str2d(val[1]);
        Dir   = str2ul(val[2]);
        Phase = str2ul(val[3]);
    }
    else if (val[0] == "PB"){   Type = IT::PB;
        Value = str2d(val[1]);
        Dir   = str2ul(val[2]);
        Phase = str2ul(val[3]);
    }
    else if (val[0] == "VB"){   Type = IT::VB;
        Value = str2d(val[1]);
        Dir   = str2ul(val[2]);
        Phase = str2ul(val[3]);
    }
    else if (val[0] == "AB")    Type = IT::AB;
    else if (val[0] == "BB")    Type = IT::BB;
    else if (val[0] == "SD"){   Type = IT::SD;
        Value = str2d(val[1]);
    }
    else if (val[0] == "HF"){   Type = IT::HF;
        Phase = str2ul(val[3]);
    }
    else throw invalid_argument("FilePrm::IDXVAL::IDXVAL() - unrecognized"
                                " index descriptor type");
}

const set<size_t> &FilePrm::getIdx(const set<IDXTYPE> &types) const
{
    size_t typesum = accumulate(types.begin(), types.end(),
                                static_cast<size_t>(0));
    try {
        return idxByTypComb.at(typesum);
    }
    catch (exception& e) {
        set<size_t> idxes;
        for (const auto& type : types) {
            for (const auto& idx : idxByType[type])
                idxes.insert(idx);
        }
        idxByTypComb[typesum] = idxes;
    }
    return idxByTypComb.at(typesum);
}

} // namespace IO

} // namespace DRL
