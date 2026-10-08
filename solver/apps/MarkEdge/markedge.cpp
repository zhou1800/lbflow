/*! \file MarkEdge/markedge.cpp
 * \author Yu Zhou (zhouyu@live.com)
 */

#include "fileprm.h"
#include "logger.h"
#include "markedge.h"

using namespace DRL;
using namespace DRL::PIPT;
using namespace DRL::HLP;
using namespace std;

/*!
 * \brief Mark Interface of solid and other with value of boundary and Save.
 *
 * Excute 2 steps:
 *  1. Load raw file.
 *  2. Mark edge based on given solid and fill.
 *  3. Save to file.
 *
 * \param prm - parameter file.
 */
template <typename T>
void MarkSave(const FilePrm& prm);

int main(int argc, char* argv[])
{
    /* Prepare thread */
    int thisId = mpi().getId();

    try {
        /* Check arguments */
        if (argc < 2){          // argv[1] is input file name
            if (thisId == 0) cerr<< argv[0] << ": you must specify the Grid"
                                               " filename." << endl;
            return 1;
        }

        /* Load input file */
        FilePrm input(argv[1]);

        /* Mark and Save */
        if (input.getFileType() == "uint8")
            MarkSave<uint8_t>(input);
        if (input.getFileType() == "uint16")
            MarkSave<uint16_t>(input);
        if (input.getFileType() == "uint32")
            MarkSave<uint32_t>(input);

    } catch (exception& e){
        cout << thisId << " Failed to excute: " << e.what() << endl;
        cout.flush();
        return -1;
    }
    /* finish */
    //cout << thisId << " Done!";
    return 0;
}

template <typename T>
void MarkSave(const FilePrm &prm) {
    Logger logger(cout);
    logger.Reset(true);

    logger.Log("Initializing...", false);
    const T SD = static_cast<T>(*(prm.getIdx({IT::SD}).begin()));
    const T BD = static_cast<T>(*(prm.getIdx({IT::BB}).begin()));
    MarkEdge marker(prm.getD(), prm.getQ());
    logger.Log("done!\n");

    logger.Log("Loading file...", false);
    valarray<T> buff = marker.LoadRaw<T>(prm.getGFile(), BD);
    logger.Log("Loaded!\n");

    logger.Log("Marking edge...", false);
    marker.Mark(&buff, SD, BD);
    logger.Log("done!\n");

    logger.Log("Saving...", false);
    string savedName = marker.SaveFile(&buff, prm.getOutPath());
    logger.Log("done!\n");

    stringstream ss;
    ss << "File saved: " << savedName << "\n"
       << "S/V = " << marker.getSV() << endl;
    logger.Log(ss.str(), false);
}
