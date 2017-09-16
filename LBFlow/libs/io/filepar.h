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
    SC          /*!< Pseudo model (shan-chen). */
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
    NA  = 0x0,  /*!< Undefined. */
    VD  = 0x1,  /*!< Void. */
    GV  = 0x2,  /*!< Void with gravity force. */
    PB  = 0x4,  /*!< Pressure boundary. */
    VB  = 0x8,  /*!< Velocity boundary. */
    SD  = 0x10, /*!< Solid. */
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
 */
typedef struct{
    char gridFile[MAX_STRING_LEN];   /*!< The gridFile. */
    int size[MAX_DIMENSION];         /*!< Grid side length. */
    int nDimension;                  /*!< Number of dimensions. */
    int nDirection;                  /*!< Number of discrete directions. */
    DECIMAL latticeLen;              /*!< Length of a unit Lattice. */
    int startTime;                   /*!< Start time. */
    int endTime;                     /*!< End time. */
    int writeInterval;               /*!< Interval for writing result file. */
    int reportInterval;              /*!< Interval for report output. */
    int restarInterval;              /*!< Interval for writing restart file. */
    DECIMAL maxError;                /*!< Terminate condition. */
    char outputPath[MAX_STRING_LEN]; /*!< Output file path name perfix. */
    int reportDirection;             /*!< Direction of report. */
    int reportSlice[NUMBER_OF_SLICE];/*!< Report slice range. */
    FLUID_MODEL fluidModel;          /*!< Fluid model. */
    EQUIL_MODEL equilModel;          /*!< Equilibrium function. */
    DECIMAL rexlationTime[BIN_PHASE];/*!< Rexlation time \f$ \tau \f$ */
    DECIMAL latticeSpeed;            /*!< Lattice speed. */
    DECIMAL alpha0;                  /*!< Free parameter. */
    DECIMAL density[BIN_PHASE];      /*!< Fluid density. */
    DECIMAL interfacialTension;      /*!< Interfacial tension \f $\sigma \f$. */
    DECIMAL interfaceThickness;      /*!< Interface thickness. */
    GRID_IDX gridIdx[MAX_GRAY];      /*!< Grid indexed.*/
} PARAMETER;

/*!
 * \brief LoadParameter loads the parameter file, then parse it.
 *
 * Store the parsed values in global variables.
 * \param[in] filePathName - parameter file path name.
 * \param[out] par - The loaded parameter structure.
 * \return NOERR on success.
 * \return ERROR on failure.
 */
int LoadParameter(const char *filePathName, PARAMETER* par);

#ifdef __cplusplus
}
#endif

#endif // FILEPAR_H
