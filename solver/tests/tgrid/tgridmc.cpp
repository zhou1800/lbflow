#include "tgrid.h"

void tGridMc(ostream& os, const GridMc& g)
{
    Debugger dbg(os);
    dbg.Add("GridMc");
    dbg.Add("AdjCi", g.AdjCi);
    dbg.Add("LCp", g.LCp);
    dbg.Add("LGs", g.LGs);
    dbg.Add("LGsWgl", g.LGsWgl);
    dbg.Add("LStart", g.LStart);
    dbg.Add("SubCp", g.getSubCp());
    dbg.Add("SubGs", g.getSubGs());
    dbg.Add("SubGsWgl", g.getSubGsWgl());
    dbg.Add("SubStart", g.getSubStart());
    dbg.Show();
    os << "-----------[GridGho.CG]-------" << endl;
    tGrid(os, g.CG);
    os << "------------------------------" << endl;
    os << "-----------[GridGho.LG]-------" << endl;
    tGrid(os, g.LG);
    os << "------------------------------" << endl;
    os << "---------[GridGho.SubG]-------" << endl;
    for(const auto& lg : g.getSubG()) tGrid(os, lg);
    os << "------------------------------" << endl;
}

void tGridMc(ostream& os, const va_idx_t& sz, size_t q, size_t ghoSz)
{
    GridMc g(sz, q, ghoSz);
    tGridMc(os, g);
}
