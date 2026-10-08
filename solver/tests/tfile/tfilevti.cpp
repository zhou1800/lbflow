#include "tfile.h"

void tFileVti(ostream& os, const FilePrm &prm)
{
    if (prm.getFileType() != "uint8"){
        os << "Test FileVti with 8-bit raw data." << endl;
        return;
    }

    FileRaw<uint8_t> raw(prm.getGFile());
    GridMc g(prm.getD(), prm.getQ());
    valarray<uint8_t> buff = raw.Load(g);
    valarray<double> mass(buff.size());
    valarray<double> velX(buff.size());
    valarray<double> velY(buff.size());
    valarray<double> velZ(buff.size());
    for (size_t i = 0; i < buff.size(); ++i){
        mass[i] = static_cast<double>(buff[i]) * 2.0;
        velX[i] = static_cast<double>(buff[i]) * 5.0;
        velY[i] = static_cast<double>(buff[i]) * 6.0;
        velZ[i] = static_cast<double>(buff[i]) * 7.0;
    }

    FileVti vti(prm.getOutPath());
    vti.AddScalar("mass", &mass);
    vti.AddVector("vel", &velX);
    vti.AddVector("vel", &velY);
    vti.AddVector("vel", &velZ);

    vti.Init(g.S.D, g.N, g.Sz, g.Stride, g.LG.Sz, g.LStart);
    vti.Save(9);
}
