/*! \file fn.h
 * \author Yu Zhou
 */
#ifndef FN_H
#define FN_H

#include "defs.h"
#include <list>

namespace DRL {

namespace ALG {

/*!
 * \brief Get the prime numbers of a given number.
 *
 * \param val - number to be decomposed.
 * \return list of prime numbers with big in the front.
 */
std::list<size_t> Prime(size_t val);

/// \name Matrix manipulation.
///@{
/*!
 * \brief Calculate the product of a matrix and a vector.
 *
 * \param mat - the matrix.
 * \param vec - the vector.
 * \param m - number of rows in the matrix.
 * \param n - number of columns in the matrix, number of rows in the vector.
 * \return a valarray of vector.
 */
va_val_t MatMulVec(const va_val_t& mat, const va_val_t& vec,
                    size_t m, size_t n);

/*!
 * \brief Calculate the product of a matrix and a diagonal matrix.
 *
 * \param mat - the matrix.
 * \param diag - the diagonal elements in the diagonal matrix.
 * \param m - number of cloumns and rows in the matrix.
 * \return a valarray of matrix.
 */
va_val_t MatMulDiag(const va_val_t& mat, const va_val_t& diag, size_t m);

/*!
 * \brief Calculate the product of two matrices.
 *
 * \param mat1 - the matrix 1.
 * \param mat2 - the matrix 2.
 * \param m - number of rows in matrix 1.
 * \param p - number of columns in matrix 1, number of rows in matrix 2.
 * \param n - number of columns in matrix 2.
 * \return a valarray of matrix.
 */
va_val_t MatMulMat(const va_val_t& mat1, const va_val_t& mat2,
                    size_t m, size_t p, size_t n);

/*!
 * \brief Claculate the transpose of a matrix.
 *
 * \param mat - the matrix.
 * \param m - number of rows in the input matrix.
 * \param n - number of columns in the input matrix.
 * \return a valarray of matrix.
 */
va_val_t MatTranspose(const va_val_t& mat, size_t m, size_t n);

/*!
 * \brief Calculate the inverse of a matrix with orthogonal columns.
 *
 * \param mat - the matrix.
 * \param m - number of rows in the matrix.
 * \param n - number of colums in the matrix.
 * \return a valarray of inversed matrix.
 */
va_val_t InvOrthoCol(const va_val_t& mat, size_t m, size_t n);

/*!
 * \brief Calculate the inverse of a matrix with orthogonal rows.
 *
 * \param mat - the matrix.
 * \param m - number of rows in the matrix.
 * \param n - number of colums in the matrix.
 * \return a valarray of inversed matrix.
 */
va_val_t InvOrthoRow(const va_val_t& mat, size_t m, size_t n);

///@}

/// \name Geometry functions.
///@{
/*!
 * \brief Generate a 2D rotate matrix.
 *
 * \param sideO - opposite side length.
 * \param sideA - adjacent side length.
 * \return a valarray of the rotate matrix.
 */
va_val_t RotateT2D(double sideO, double sideA);

/*!
 * \brief Generate a 3D rotate matrix around X axis.
 *
 * \param sideO - opposite side length.
 * \param sideA - adjacent side length.
 * \return a valarray of the rotate matrix.
 */
va_val_t RotateTx3D(double sideO, double sideA);

/*!
 * \brief Generate a 3D rotate matrix around Y axis.
 *
 * \param sideO - opposite side length.
 * \param sideA - adjacent side length.
 * \return a valarray of the rotate matrix.
 */
va_val_t RotateTy3D(double sideO, double sideA);

/*!
 * \brief Generate a 3D rotate matrix around Z axis.
 *
 * \param sideO - opposite side length.
 * \param sideA - adjacent side length.
 * \return a valarray of the rotate matrix.
 */
va_val_t RotateTz3D(double sideO, double sideA);

///@}

} // namespace ALG

} // namespace DRL

#endif // FN_H
