/*! \file fileprm.cpp
 * \author Yu Zhou
 */
#include "fileprm.h"
#include "thread.h"
#include <fstream>
#include <sstream>
#include <cstring>
#include <numeric>

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
