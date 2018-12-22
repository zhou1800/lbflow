/*! \file config.h
 * \brief Global configurations.
 */

#ifndef CONFIG_H
#define CONFIG_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#define DECIMAL float       /*!< Decimal value type. */
#define RAW unsigned char   /*!< Raw file type. */

/*!
 * \brief The MODEL_CONSTS enum defines constants for the simulator.
 */
enum MODEL_CONSTS{
    MAX_DIMENSION = 3,      /*!< Maximum number of dimensions. */
    MAX_GRAYSCALE = 256,    /*!< Grid file grayscale range: 0 ~ 255. */
    ROI_SLICE = 6,          /*!< Number of ROI slice positions (x+, x-, y+, y-, z+, z-). */
    SINGLE_PHASE = 1,       /*!< Single phase. */
    BIN_PHASE = 2,          /*!< Binary phase. */
};

enum PROGRAM_STATUS{
    SUCCESS = 0,            /*!< Return value on successful function call. */
    ERROR = -1,             /*!< Return value on function call failure. */
};

/*!
 * If expr is TRUE, print an error message and goto end.
 * \param[in] expr - expression.
 * \param[in] ... - error message.
 */
#define TO_END(expr, ...)\
    if ((expr)) { fprintf(stderr, __VA_ARGS__); goto end; }
/*!
 * If expr is FALSE, print an error message and terminate functiona (return ERROR).
 * \param[in] expr - expression.
 * \param[in] ... - error message.
 */
#define ASSERT(expr, ...)\
    if (!(expr)) { fprintf(stderr, __VA_ARGS__); return ERROR; }

/*!
 * If pointer is not NULL, free the pointer and set it to NULL.
*/
#define FREE(pointer)\
    { free(pointer); pointer=NULL; }

#ifdef __cplusplus
}
#endif

#endif // CONFIG_H
