#include "tfile.h"
#include <stdexcept>
#include <sstream>
using std::exception;

void tFilePrm(ostream& os, const string &fname)
{
    FilePrm f(fname);

    Debugger dbg(os);
    dbg.Add("FilePrm");
    try                  {dbg.Add("Alpha0", f.getAlpha0());}
    catch (exception& e) {os << e.what() << endl; }
    try                  {dbg.Add("cAngle", f.getCAngle());}
    catch (exception& e) {os << e.what() << endl; }
    try                  {dbg.Add("Components", f.getComponents());}
    catch (exception& e) {os << e.what() << endl; }
    try                  {dbg.Add("Density", f.getDensity());}
    catch (exception& e) {os << e.what() << endl; }
    try                  {dbg.Add("Dimension", f.getD());}
    catch (exception& e) {os << e.what() << endl; }
    try                  {dbg.Add("Directions", f.getQ());}
    catch (exception& e) {os << e.what() << endl; }
    try                  {dbg.Add("FileType", f.getFileType());}
    catch (exception& e) {os << e.what() << endl; }
    try                  {dbg.Add("FluidModel", f.getFluidModel());}
    catch (exception& e) {os << e.what() << endl; }
    try                  {dbg.Add("GridFile", f.getGFile());}
    catch (exception& e) {os << e.what() << endl; }
    try                  {dbg.Add("IThickness", f.getIThickness());}
    catch (exception& e) {os << e.what() << endl; }
    try                  {dbg.Add("LatticeSpeed", f.getLSpeed());}
    catch (exception& e) {os << e.what() << endl; }
    try                  {dbg.Add("LogPath", f.getLogPath());}
    catch (exception& e) {os << e.what() << endl; }
    try                  {dbg.Add("OutPath", f.getOutPath());}
    catch (exception& e) {os << e.what() << endl; }
    try                  {dbg.Add("Phase", f.getPhase());}
    catch (exception& e) {os << e.what() << endl; }
    try                  {dbg.Add("RexlationTime", f.getRTime());}
    catch (exception& e) {os << e.what() << endl; }
    try                  {dbg.Add("Tension", f.getSTension());}
    catch (exception& e) {os << e.what() << endl; }
    try                  {dbg.Add("TEnd", f.getTEnd());}
    catch (exception& e) {os << e.what() << endl; }
    try                  {dbg.Add("TRst", f.getTRst());}
    catch (exception& e) {os << e.what() << endl; }
    try                  {dbg.Add("TStart", f.getTStart());}
    catch (exception& e) {os << e.what() << endl; }
    try                  {dbg.Add("TStep", f.getTStep());}
    catch (exception& e) {os << e.what() << endl; }

    vector<string> names;
    vector<string> values;
    try {
        std::stringstream ss;
        for (auto a :
             f.getIdx({IT::PB, IT::VD, IT::GV, IT::BB, IT::VB, IT::SD})){
            ss << "Index" << a << " ";
            names.push_back(ss.str());
            ss.str("");
            ss << f.getIVal(a).Type << " "
               << f.getIVal(a).Value << " "
               << f.getIVal(a).Dir  << " "
               << f.getIVal(a).Phase;
            values.push_back(ss.str());
            ss.str("");
        }
    }
    catch (exception& e) {
        os << e.what() << endl;
    }
    for (size_t i = 0; i < names.size(); ++i)
        dbg.Add(names[i], values[i]);
    dbg.Show();
}
