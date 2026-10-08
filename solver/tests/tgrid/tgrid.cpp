#include "tgrid.h"

void tGrid(ostream& os, const va_idx_t& sz, size_t q)
{
    Grid g(sz, q);
    tGrid(os, g);
}

void tGrid(ostream& os, const Grid& g)
{
    Debugger dbg(os);
    dbg.Add("Grid");
    dbg.Add("N", g.N);
    dbg.Add("Stride", g.Stride);
    dbg.Add("Sz", g.Sz);
    dbg.Add("AdjI", g.getAdjI());
    dbg.Show();
    os << "---------------[Grid.S]-------" << endl;
    tGridSpace(os, g.S);
    os << "------------------------------" << endl;
}
