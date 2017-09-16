/*! \file filepar.h
 * \brief Read and parse the parameter file.
 * \author Yu Zhou (zhouyu@live.com)
 */
#ifndef FILEPAR_H
#define FILEPAR_H

#include "config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! Conventions of parameter file. */
enum PAR_FILE{
    MAX_GRAY = 256,         /*!< Grid file grayscale range: 0 ~ 255. */
    MAX_STRING_LEN = 256,   /*!< Maximum number of characters in a string. */
    MAX_KEY_LEN = 32,       /*!< Maximum number of characters in a keyword. */
    MAX_FILE_LEN = 8192,    /*!< Maximum length of parameter file. */
};

/*! Fluid model. */
typedef enum {
    SCSP,       /*!< Single phase single components model. */
    COLOR,      /*!< Binary color fluid model. */
} FLUID_MODEL;

/*! Equilibrium function. */
typedef enum {
    SRT,        /*!< Single relaxation time. */
    MRT         /*!< Mutiple relaxation time. */
} EQUIL_MODEL;

/*!
 * \brief GRIDTYPE defines the type of grid.
 *
 * The value of grid types are assigned by \f$ 2^n \f$ to get a unique
 *  sum value of any type combination.
 */
typedef enum {
    NA  = (0),       /*!< Undefined. */
    VD  = (1 << 0),  /*!< Void. */
    GV  = (1 << 1),  /*!< Void with gravity force. */
    PB  = (1 << 2),  /*!< Pressure boundary. */
    VB  = (1 << 3),  /*!< Velocity boundary. */
    SD  = (1 << 4),  /*!< Solid. */
} GRID_TYPE;

/*!
 * \brief The GRID_IDX struct is the index of grid.
 *
 * This struct stores the attributes corresponding to each grayscale value,
 *  which is defined in the parameter file.
 *
 * \note (32bit, MinGW) DECIMAL = double, size = 24; DECIMAL = float, size = 16.
 */
typedef struct {
    GRID_TYPE type;         /*!< Grid type. */
    DECIMAL value;          /*!< Value of the force. */
    int dir;                /*!< Force direction. */
    int phase;              /*!< Phase index. */
} GRID_IDX;

/*!
 * \brief The PARAMETER struct is an image of the parameter file.
 *
 * \note (32bit, MinGW) DECIMAL = double, size = 6552; DECIMAL = float, size =
 *  4452.
 * \note Always initialize with DEFAULT_PARAMETER.
 */
typedef struct{
    char gridFile[MAX_STRING_LEN];  /*!< The gridFile. */
    int size[MAX_DIMENSION];        /*!< Grid side length. */
    int nDimension;                 /*!< Number of dimensions. */
    int nDirection;                 /*!< Number of discrete directions. */
    DECIMAL latticeLen;             /*!< Length of a unit Lattice. */
    int startTime;                  /*!< Start time. */
    int endTime;                    /*!< End time. */
    int writeIntvl;                 /*!< Interval for writing result file. */
    int rptIntvl;                   /*!< Interval for report output. */
    int rstIntvl;                   /*!< Interval for writing restart file. */
    DECIMAL maxError;               /*!< Terminate condition. */
    char outPath[MAX_STRING_LEN];   /*!< Output file path name perfix. */
    int rptDir;                     /*!< Direction of report. */
    int rptSlice[NUMBER_OF_SLICE];  /*!< Report slice range. */
    FLUID_MODEL fluidModel;         /*!< Fluid model. */
    EQUIL_MODEL equilModel;         /*!< Equilibrium function. */
    DECIMAL rlxTime[BIN_PHASE];     /*!< Relaxation time \f$ \tau \f$ */
    DECIMAL latticeSpd;             /*!< Lattice speed. */
    DECIMAL density[BIN_PHASE];     /*!< Fluid density. */
    DECIMAL alpha0;                 /*!< Free parameter. */
    DECIMAL ift;                    /*!< Interfacial tension \f $\sigma \f$. */
    DECIMAL intfThk;                /*!< Interface thickness. */
    GRID_IDX index[MAX_GRAY];     /*!< Grid index.*/
} PARAMETER;

/*!
 * \brief LoadParameter loads the parameter file, then parse it.
 *
 * Store the parsed values in global variables.
 * \param[in] filePathName - parameter file path name.
 * \param[out] par - The loaded parameter structure.
 * \return NOERR on success. exit(ERROR) on failure.
 */
int LoadParameter(const char *filePathName, PARAMETER* par);

#ifdef __cplusplus
}
#endif

#endif // FILEPAR_H
