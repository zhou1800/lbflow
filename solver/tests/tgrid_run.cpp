#include "tgrid.h"
#include "thread.h"
#include <sstream>
#include <fstream>

using namespace std;
using namespace DRL;

void TestGridSpace(ostream& os);
void TestGrid(ostream& os);
void TestGridGho(ostream& os);
void TestGridMc(ostream& os);
void TestROI(ostream& os);

int main()
{
    try {
        ostringstream os;
        if (mpi().getId() == 0){
            os << "TestGridSpace.txt";
            ofstream fs(os.str());
            TestGridSpace(fs);
            os.str("");
            fs.close();
            os << "TestGrid.txt";
            fs.open(os.str());
            TestGrid(fs);
            os.str("");
            fs.close();
            os << "TestGridGho.txt";
            fs.open(os.str());
            TestGridGho(fs);
            os.str("");
            fs.close();
            os << "TestROI.txt";
            fs.open(os.str());
            TestROI(fs);
            os.str("");
            fs.close();
        }

        os << "TestGridMc" << mpi().getId() << ".txt";
        ofstream fs(os.str());
        TestGridMc(fs);
        fs.close();
    }
    catch (exception &e)
    {
        cout << mpi().getId() << " Failed to excute: " << e.what() << endl;
        cout.flush();
        return -1;
    }

}

void TestGridSpace(ostream &os)
{
    tGridSpace(os, 1, 3);
    tGridSpace(os, 2, 5);
    tGridSpace(os, 2, 9);
    tGridSpace(os, 3, 19);
    tGridSpace(os, 3, 27);
}

void TestGrid(ostream &os)
{
    tGrid(os, {2}, 3);
    tGrid(os, {2,3}, 9);
    tGrid(os, {2,3,4}, 19);
}

void TestGridGho(ostream &os)
{
    tGridGho(os, {2}, 3, 5);
    tGridGho(os, {2,3}, 9, 2);
    tGridGho(os, {2,3,4}, 19, 1);
}

void TestGridMc(ostream& os)
{
    tGridMc(os, {10}, 3, 5);
    tGridMc(os, {3,4}, 9, 2);
    tGridMc(os, {2,3,4}, 19, 1);
}

void TestROI(ostream& os)
{
    GridGho g({4,2,3}, 19);
    va_u8_t raw(24);
    for (size_t i = 0; i < 24; ++i)
        raw[i] = i%3;
    tROI(os, g, raw);
}
