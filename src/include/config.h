/*! \file config.h
 * \brief Global configurations.
 */

#ifndef CONFIG_H
#define CONFIG_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#define DECIMAL float       /*!< Precision of the decimal value. */
#define RAW unsigned char   /*!< Raw file type. */

/*!
 * \brief The MODEL_CONSTS enum defines constants for the simulator.
 */
enum MODEL_CONSTS{
    MAX_DIMENSION = 3,      /*!< Maximum number of dimensions. */
    ROI_SLICE = 6,          /*!< Number of ROI slice positions (x+, x-, y+, y-, z+, z-). */
    SINGLE_PHASE = 1,       /*!< Single phase. */
    BIN_PHASE = 2,          /*!< Binary phase. */
};

enum PROGRAM_STATUS{
    SUCCESS = 0,            /*!< Return value on successful function call. */
    ERROR = -1,             /*!< Return value on function call failure. */
};

/*!
 * If expr is FALSE, print an error message and return ERROR.
 * \param[in] expr - expression.
 * \param[in] ... - error message.
 */
#define VALIDATE(expr, ...)\
    if (!(expr)) { fprintf(stderr, __VA_ARGS__); return ERROR; }

#ifdef __cplusplus
}
#endif

#endif // CONFIG_H
