#include "tgrid.h"

void tGridSpace(ostream& os, size_t d, size_t q)
{
    GridSpace gsp(d, q);
    tGridSpace(os, gsp);
}

void tGridSpace(ostream &os, const GridSpace &gsp)
{
    Debugger dbg(os);
    dbg.Add("GridSpace");
    dbg.Add("D", gsp.D);
    dbg.Add("Q", gsp.Q);
    dbg.Add("Qi", gsp.Qi);
    dbg.Add("Qp", gsp.Qp);
    dbg.Add("Qop", gsp.Qop);
    dbg.Add("Qhalf", gsp.Qhalf);
    dbg.Add("Qadj", gsp.Qadj);
    dbg.Add("Qrel", gsp.Qrel);
    dbg.Add("Qnorm", gsp.Qprep);
    dbg.Add("Qaxis", gsp.Qaxis);
    dbg.Add("Qsign", gsp.Qsign);
    dbg.Show();
}
