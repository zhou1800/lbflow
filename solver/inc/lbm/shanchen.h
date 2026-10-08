/*! \file shanchen.h
 * \author Yu Zhou
 */
#ifndef SHANCHEN_H
#define SHANCHEN_H

#include "simulator.h"

namespace DRL {

namespace LBM {

class ShanChen : public Simulator
{
public:
    ShanChen(const FilePrm& prm);
    virtual ~ShanChen();

    virtual void Run(std::ostream &logProc = cout) override;
};

} // namespace LBM

} // namespace DRL

#endif // SHANCHEN_H
