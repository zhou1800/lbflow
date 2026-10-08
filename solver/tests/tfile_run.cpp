#include "tfile.h"
#include "thread.h"
#include <sstream>
#include <fstream>

using namespace std;
using namespace DRL;

void TestFile(ostream& os, const string& fname);
void TestFilePrm(ostream& os, const string &fname);
void TestFileRaw(ostream& os, const string& fname);
void TestFileVti(ostream& os, const string& fname);

int main(int argc, char* argv[])
{
    try {
        /* check program arguments */
        if (argc < 2){          // argv[1] is input file name
            if (mpi().getId() == 0)
                cerr << argv[0]
                        << ": you must specify an input filename." << endl;
            return 1;
        }

        if (mpi().getId() == 0)
            TestFile(cout , argv[1]);

        ostringstream os;
        os << "TestFilePrm" << mpi().getId() << ".txt";
        ofstream fs(os.str());

        TestFilePrm(fs, argv[1]);
        TestFileRaw(fs, argv[1]);
        TestFileVti(fs, argv[1]);
        fs.close();
    }
    catch (exception &e)
    {
        cout << mpi().getId() << " Failed to excute: " << e.what() << endl;
        cout.flush();
        return -1;
    }

}

void TestFile(ostream& os, const string &fname)
{
    tFile(os, fname);
    tFile(os, "[fake file]");
}

void TestFilePrm(ostream& os, const string &fname)
{
    tFilePrm(os, fname);
    try {
        tFilePrm(os, "[fake file]");
    }
    catch (exception& e) {
        os << "[fake file] " << e.what() << endl;
    }
}

void TestFileRaw(ostream& os, const string &fname)
{
    FilePrm prm(fname);
    if (prm.getFileType() == "uint8")
        tFileRaw<uint8_t>(os, prm);
    if (prm.getFileType() == "uint16")
        tFileRaw<uint16_t>(os, prm);
    if (prm.getFileType() == "uint32")
        tFileRaw<uint32_t>(os, prm);
}

void TestFileVti(ostream& os, const string &fname)
{
    FilePrm prm(fname);
    tFileVti(os, prm);
}
