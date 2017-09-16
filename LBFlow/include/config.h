/*! \file config.h
 * \brief Global configurations.
 * \author Yu Zhou
 */

#ifndef CONFIG_H
#define CONFIG_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#define DECIMAL float       /*!< Precision of the decimal value. */

/*!
 * \brief The MODEL_CONSTS enum defines constants for the simulator.
 */
enum MODEL_CONSTS{
    MAX_DIMENSION = 3,      /*!< Maximum number of dimensions. */
    NUMBER_OF_SLICE = 2,    /*!< Number of report slices. */
    BIN_PHASE = 2,          /*!< Number of binary phase. */
};

enum PROGRAM_STATUS{
    NOERR = 0,          /*!< Return value on successful function call. */
    ERROR = -1,         /*!< Return value on function call failure. */
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
