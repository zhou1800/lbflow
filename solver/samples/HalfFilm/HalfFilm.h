    #grid
    {
        gridFile     "23x12y40z.raw";
        dimension    40 12 23;
        directions   19;
        fileType     "uint8";
        latticeLen   0.000001;
    }
    #control
    {
        startTime       0;
        endTime         10000;
        writeInterval   1000;
        reportInterval  1000;
        restartInterval 10000;
        stopError       0.000001;
        outPath         "./COLOR/HalfFilm";
        logPath         ".";
        reportDir       5;
        reportSlice     10 38;
    }
    #scsp
    {
        //fluidModel          "scsp";
        //equilModel          "SRT";
        //rexlationTime       1;
        //latticeSpeed        1;
        //density             1;
    }
    #colorfluid
    {
        fluidModel          "color";
        equilModel          "MRT";
        rexlationTime       0.57    1.215416667;
        latticeSpeed        1;
        alpha0              0.33333333333333333333;
        density             1       0.831683168;
        contactAngle        180;
        surfaceTension      0.085920139;
        interfaceThickness  0.7;
    }
    #SC
    {
        //fluidModel          "shanchen"
        //components  "oil"   "water";
    }
    #Index
    {
      //        0       1             2       3
      //index   type    value         dir     phase
        0       "VD"    "N/A"         "N/A"   "0"
        25      "HF"    "N/A"         "N/A"   "0"
        50      "VB"    "0.0001"      "5"     "0"
        100     "VD"    "N/A"         "N/A"   "1"
        125     "HF"    "N/A"         "N/A"   "1"
        150     "VB"    "0.0001"      "5"     "1"
        200     "PB"    "1"           "6"     "0"
        255     "SD"    "N/A"         "N/A"   "N/A"
    }
