#include "tgrid.h"

void tGridGho(ostream& os, const GridGho& g)
{
    Debugger dbg(os);
    dbg.Add("GridGho");
    dbg.Add("i2fGs", g.getI2fGs());
    dbg.Add("i2fMap", g.getI2fMap());
    dbg.Add("GhoSz", g.GhoSz);
    dbg.Show();
    os << "-----------[GridGho.FG]-------" << endl;
    tGrid(os, g.FG);
    os << "------------------------------" << endl;
}

void tGridGho(ostream& os, const va_idx_t& sz, size_t q, size_t ghoSz)
{
    GridGho g(sz, q, ghoSz);
    tGridGho(os, g);
}
