/*! \file fn.cpp
 * \author Yu Zhou
 */
#include "fn.h"
#include <stdexcept>
#include <numeric>

namespace DRL {

namespace ALG {

using namespace std;

list<size_t> Prime(size_t val)
{
    list<size_t> primes;
    size_t i = 2;
    while (val > 1) {
        if (val % i == 0) {
            primes.push_front(i);
            val /= i;
        }
        else
            i++;
    }
    return move(primes);
}

va_val_t MatMulVec(const va_val_t &mat, const va_val_t &vec, size_t m, size_t n)
{
    va_val_t r(0.0, m);        // result
    auto imat = begin(mat);
    for (size_t i = 0; i < m; ++i){
        r[i] = inner_product(begin(vec), end(vec), imat,
                             static_cast<double>(0.0));
        imat += n;
    }
    return move(r);
}

va_val_t MatMulDiag(const va_val_t &mat, const va_val_t &diag, size_t m)
{
    va_val_t r(0.0, m * m);        // result
    auto imat = begin(mat);
    auto ir = begin(r);
    for (size_t i = 0; i < m; ++i){
        for (size_t j = 0; j < m; ++j)
            *ir++ = *imat++ * diag[j];
    }
    return move(r);
}

va_val_t MatMulMat(const va_val_t &mat1, const va_val_t &mat2,
                    size_t m, size_t p, size_t n)
{
    va_val_t r(0.0, m * n);        // result
    va_val_t m2 = MatTranspose(mat2, p, n);    // transposed mat2
    auto i1 = begin(mat1);
    auto ir = begin(r);
    for (size_t i = 0; i < m; ++i){
        auto i2 = begin(m2);
        for (size_t j = 0; j < n; ++j){
            *ir++ = inner_product(i1, i1 + p, i2, 0.0);
            i2 += p;
        }
        i1 += p;
    }
    return move(r);
}

va_val_t MatTranspose(const va_val_t &mat, size_t m, size_t n)
{
    va_val_t r(0.0, m * n);        // result
    auto ir = begin(r);
    for (size_t i = 0; i < m; ++i){
        auto im = begin(mat) + i;
        for (size_t j = 0; j < n; ++j){
            *ir++ = *im;
            im += n;
        }
    }
    return move(r);
}

va_val_t InvOrthoCol(const va_val_t &mat, size_t m, size_t n)
{
    va_val_t r(0.0, mat.size());   // result
    for (size_t j = 0; j < n; ++j){     // for columns
        double s = 0;
        for (size_t i = 0; i < m; ++i) {
            size_t idx = n*i + j;
            s += mat[idx] * mat[idx];
        }
        s = 1/s;
        size_t iS = j*m;
        for (size_t i = 0; i < m; ++i)
            r[iS + i] = mat[n*i + j] * s;
    }
    return move(r);
}

va_val_t InvOrthoRow(const va_val_t &mat, size_t m, size_t n)
{
    va_val_t r(0.0, mat.size());   // result
    auto imat = begin(mat);
    for (size_t j = 0; j < m; ++j){     // for rows of origional
        auto imatEnd = imat+n;
        double s = inner_product(imat, imatEnd, imat, 0.0);
        imat = imatEnd;
        s = 1/s;
        size_t iS = j*n;
        for (size_t i = 0; i < n; ++i)  // cor columns of origional
            r[m*i + j] = mat[iS + i] * s;
    }
    return move(r);
}

va_val_t RotateT2D(double sideO, double sideA)
{
    double c = 1;
    double s = 0;
    double hypotenuse = sqrt(sideO * sideO + sideA * sideA);
    if (hypotenuse != 0){
        c = sideA / hypotenuse;
        s = sideO / hypotenuse;
    }
    return {c, -s, s, c};
}

va_val_t RotateTx3D(double sideO, double sideA)
{
    double c = 1;
    double s = 0;
    double hypotenuse = sqrt(sideO * sideO + sideA * sideA);
    if (hypotenuse != 0){
        c = sideA / hypotenuse;
        s = sideO / hypotenuse;
    }
    return {c, -s, 0, s, c, 0, 0, 0, 1};
}

va_val_t RotateTy3D(double sideO, double sideA)
{
    double c = 1;
    double s = 0;
    double hypotenuse = sqrt(sideO * sideO + sideA * sideA);
    if (hypotenuse != 0){
        c = sideA / hypotenuse;
        s = sideO / hypotenuse;
    }
    return {c, 0, s, 0, 1, 0, -s, 0, c};
}

va_val_t RotateTz3D(double sideO, double sideA)
{
    double c = 1;
    double s = 0;
    double hypotenuse = sqrt(sideO * sideO + sideA * sideA);
    if (hypotenuse != 0){
        c = sideA / hypotenuse;
        s = sideO / hypotenuse;
    }
    return {1, 0, 0, 0, c, -s, 0, s, c};
}

} // namespace ALG

} // namespace DRL

