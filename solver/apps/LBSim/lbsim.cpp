/*! \file lbsim.cpp
 * \author Yu Zhou (zhouyu@live.com)
 */

#include "logger.h"
#include "colorfluid.h"
#include "scsp.h"
#include "shanchen.h"

using namespace std;
using namespace DRL;
using namespace DRL::IO;
using namespace DRL::LBM;
using namespace DRL::HLP;

int main(int argc, char* argv[])
{
    /* prepare thread */
    int thisId = mpi().getId();

    try {
        /* check program arguments */
        if (argc < 2){          // argv[1] is input file name
            if (thisId == 0)
                cerr<< argv[0]
                        << ": you must specify an input filename." << endl;
            return 1;
        }

        /* load input file */

        Logger logger(cout);

        logger.Log("Model constructing...", false);
        FilePrm input(argv[1]);

        Simulator* simulator;

        // init fluid model
        if (input.getFluidModel() == "color")
            simulator = new ColorFluid(input);
        else if (input.getFluidModel() == "scsp") {
            simulator = new SCSP(input);
        }
        else if (input.getFluidModel() == "shanchen") {
            simulator = new ShanChen(input);
            throw invalid_argument("Main(): Fluid model shanchen currently not "
                                   "supported.");
        }
        else
            throw invalid_argument("Main(): Unsupport fluid model (only color "
                                   "model and single components single phase "
                                   "model are supported at this time)");

        mpi().Barrier();
        logger.Log("Done!\n");

        logger.Reset();
        logger.Log("Simulation started!\n",false);

        simulator->Run(cout);

        mpi().Barrier();
        logger.Log("Simulation finished!\nTotal ");

    } catch (exception& e){
        cout << thisId << " Failed to excute: " << e.what() << endl;
        cout.flush();
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        return EXIT_FAILURE;
    }
    /* finish */
    return 0;
}
