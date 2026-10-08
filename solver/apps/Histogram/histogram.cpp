/*! \file Histogram/histogram.cpp
 * \author Yu Zhou
 */

#include "fileprm.h"
#include "logger.h"
#include "histogram.h"
#include <sstream>

using namespace DRL;
using namespace DRL::HLP;
using namespace DRL::PIPT;
using namespace std;

/*!
 * \brief Calculate histogram and save to file.
 *
 * Excute 3 steps:
 *  1. Load raw file.
 *  2. Calculate histogram.
 *  4. Save histogram to file.
 *
 * \param prm - parameter file.
 */
template <typename T>
void HistSave(const FilePrm& prm);

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

        /* Call Calculate and Save */
        if (input.getFileType() == "uint8")
            HistSave<uint8_t>(input);
        else if (input.getFileType() == "uint16")
            HistSave<uint16_t>(input);
        else if (input.getFileType() == "uint32")
            HistSave<uint32_t>(input);

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
void HistSave(const FilePrm& prm)
{
    Logger logger(cout);
    logger.Reset(true);

    logger.Log("Initializing...", false);
    GridMc g(prm.getD(), prm.getQ());
    Histogram hist(g);
    logger.Log("done!\t");

    logger.Log("Loading file...", false);
    valarray<T> buff = hist.LoadRaw<T>(prm.getGFile());
    logger.Log("Loaded!\n");

    logger.Log("Calculating histogram...", false);
    hist.CalHist(buff);
    logger.Log("Done!\n");

    string savedName = hist.SaveHistFile(prm.getLogPath(), prm.getGFile());
    if (mpi().getId() == 0) {
        std::stringstream ss;
        ss << "Histogram file saved: " << savedName << "\n";
        logger.Log(ss.str(), false);
    }
}
