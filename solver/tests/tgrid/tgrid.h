#ifndef TGRID_H
#define TGRID_H

#include "debugger.h"
#include "gridspace.h"
#include "grid.h"
#include "gridgho.h"
#include "gridmc.h"

using namespace DRL;
using namespace DRL::TEST;
using namespace DRL::GRID;

void tGridSpace(ostream &os, const GridSpace& gsp);
void tGridSpace(ostream &os, size_t d, size_t q);
void tGrid(ostream& os, const Grid& g);
void tGrid(ostream& os, const va_idx_t& sz, size_t q);
void tGridGho(ostream& os, const GridGho& g);
void tGridGho(ostream& os, const va_idx_t& sz, size_t q, size_t ghoSz);
void tGridMc(ostream& os, const GridMc& g);
void tGridMc(ostream& os, const va_idx_t& sz, size_t q, size_t ghoSz);
void tROI(ostream& os, const GridGho& g, const va_u8_t& raw);

#endif // TGRID_H
