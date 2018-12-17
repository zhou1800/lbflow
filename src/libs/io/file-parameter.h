/*! \file file-parameter.h
 * \brief Read and parse the parameter file.
 */
#ifndef FILE_PARAMETER_H
#define FILE_PARAMETER_H

#include "config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! Conventions of parameter file. */
enum PARAMETE_FILE{
    MAX_GRAYSCALE = 256,    /*!< Grid file grayscale range: 0 ~ 255. */
    MAX_STRING_LEN = 256,   /*!< Maximum number of characters in a string. */
    MAX_KEY_LEN = 32,       /*!< Maximum number of characters in a keyword. */
    MAX_FILE_LEN = 8192,    /*!< Maximum length of parameter file. */
};

/*! Fluid model. */
typedef enum {
    SCSP,       /*!< Single phase single components model. */
    COLOR,      /*!< Binary color fluid model. */
} FLUID_MODEL;

/*! Relaxation model. */
typedef enum {
    SRT,        /*!< Single relaxation time. */
    MRT         /*!< Mutiple relaxation time. */
} RELAX_MODEL;

/*!
 * \brief GRID_TYPE defines the type of grid.
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
 * \brief The GRAYSCALE struct defines the property of grayscale value in the grid.
 *
 * This struct stores the attributes corresponding to each grayscale value,
 *  which is defined in the parameter file.
 */
typedef struct {
    GRID_TYPE type;         /*!< Grid type. */
    DECIMAL value;          /*!< Value of the physical quantity. */
    int direction;          /*!< Direction of the physical quantity. */
    int phase;              /*!< Phase index. */
} GRAYSCALE;

/*!
 * \brief The PARAMETER struct is an image of the parameter file.
 * \note Always initialize with DEFAULT_PARAMETER.
 */
typedef struct{
    /* Grid */
    char gridFile[MAX_STRING_LEN];      /*!< The gridFile. */
    int size[MAX_DIMENSION];            /*!< Grid side length. */
    int nDimension;                     /*!< Number of dimensions. */
    int nDirection;                     /*!< Number of discrete directions. */
    /* Control */
    int startStep;                      /*!< Start time step. */
    int endStep;                        /*!< End time step. */
    int resultIntvl;                    /*!< Interval for writing result file. */
    int restartIntvl;                   /*!< Interval for writing restart file. */
    int reportIntvl;                    /*!< Interval for report output. */
    DECIMAL maxError;                   /*!< Terminate condition. */
    /* Output */
    char outPath[MAX_STRING_LEN];       /*!< Output file dir with name perfix. */
    int reportROI[ROI_SLICE];           /*!< ROI for report. Where `slice==-1` is not reported.*/
    /* Model */
    FLUID_MODEL fluidModel;             /*!< Fluid model. */
    RELAX_MODEL relaxModel;             /*!< Relaxation model. */
    DECIMAL density[BIN_PHASE];         /*!< Fluid density. */
    DECIMAL relaxTime[BIN_PHASE];       /*!< Relaxation time \f$ \tau \f$ */
    DECIMAL latticeSpeed;               /*!< Lattice speed. */
    DECIMAL alpha0;                     /*!< Free parameter. */
    DECIMAL ift;                        /*!< Interfacial tension \f$\sigma \f$. */
    DECIMAL thickness;                  /*!< Interface thickness. */
    GRAYSCALE grayscale[MAX_GRAYSCALE]; /*!< Grayscale defination.*/
} PARAMETER;

/*!
 * \brief The LoadParameter function loads the parameter file, and then parse it.
 *
 * \param[in] pathName - parameter file path name.
 * \param[out] par - The loaded parameter structure.
 * \return SUCCESS on success. exit(ERROR) on failure.
 */
int LoadParameter(const char *pathName, PARAMETER* par);

#ifdef __cplusplus
}
#endif

#endif // FILE_PARAMETER_H
