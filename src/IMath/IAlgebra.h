#ifndef IALGEBRA_H
#define IALGEBRA_H

#include "IFunc.h"
#include "IVector.h"
#include "IMatrix.h"

#include <cmath>
#include <cstddef>
#include <algorithm>
#include <vector>
#include <type_traits>


namespace IMath
{

/* --- Global Functions --- */

//! Returns the value of 1 + 2 + ... + n = n*(n+1)/2.
template <typename T>
T GaussianSum(T n)
{
    static_assert(std::is_integral<T>::value, "GaussianSum function only allows integral types");
    return n*(n + T(1))/T(2);
}

//! Returns the value of 1^2 + 2^2 + ... + n^2 = n*(n+1)*(2n+1)/6.
template <typename T>
T GaussianSumSq(T n)
{
    static_assert(std::is_integral<T>::value, "GaussianSumSq function only allows integral types");
    return n*(n + T(1))*(n*T(2) + T(1))/T(6);
}

//! Computes a normal (gaussian) distribution value for the specified (1-dimensional) position x with the specified mean and variance.
template <typename T>
T NormalDistribution(const T& x, const T& mean, const T& variance)
{
    return std::exp(-(x - mean)*(x - mean) / (variance + variance)) / std::sqrt(T(2) * T(M_PI) * variance);
}

//! Computes a normal (gaussian) distribution value for the specified (1-dimensional) position x with the specified mean and variance.
template <typename T>
T NormalDistribution(const T& x)
{
    return std::exp(-(x*x) / T(2)) / std::sqrt(T(2) * T(M_PI));
}


//! Returns the angle (in radians) between the two normalized vectors 'lhs' and 'rhs'.
template <typename VectorType, typename ScalarType = typename VectorType::ScalarType>
ScalarType AngleNorm(const VectorType& lhs, const VectorType& rhs)
{
    return std::acos(Dot<VectorType, ScalarType>(lhs, rhs));
}


//! Returns the reflected vector of the incident vector for the specified surface normal.
template <typename VectorType, typename ScalarType = typename VectorType::ScalarType>
VectorType Reflect(const VectorType& incident, const VectorType& normal)
{
    /* Compute reflection as: I - N x Dot(N, I) x 2 */
    auto v = normal;
    v *= (Dot<VectorType, ScalarType>(normal, incident) * ScalarType(-2));
    v += incident;
    return v;
}


/**
\brief Mixes the two values by their scalings.
\return Equivalent to: v0*scale0 + v1*scale1
*/
template <typename T, typename I>
T Mix(const T& v0, const T& v1, const I& scale0, const I& scale1)
{
    return v0*scale0 + v1*scale1;
}

/**
\brief Clamps the input value 'x' into the range [0, 1].
\return max{ 0, min{ x, 1 } }
*/
template <typename T>
T Saturate(const T& x)
{
    return std::max(T(0), std::min(x, T(1)));
}

/**
\brief Clamps the value 'x' into the range [minima, maxima].
\return max{ minima, min{ x, maxima } }
*/
template <typename T>
T Clamp(const T& x, const T& minima, const T& maxima)
{
    if (x <= minima)
        return minima;
    if (x >= maxima)
        return maxima;
    return x;
}

/**
\brief Returns the spherical linear interpolation between the two quaternions 'from' and 'to'.
\see QuaternionT::Slerp
*/
template <typename VectorType, typename ScalarType = typename VectorType::ScalarType>
VectorType Slerp(const VectorType& from, const VectorType& to, const ScalarType& t)
{
    ScalarType omega, cosom, sinom;
    ScalarType scale0, scale1;

    /* Calculate cosine */
    cosom = Dot<VectorType, ScalarType>(from, to);

    /* Adjust signs (if necessary) */
    if (cosom < ScalarType(0))
    {
        cosom = -cosom;
        scale1 = ScalarType(-1);
    }
    else
        scale1 = ScalarType(1);

    /* Calculate coefficients */
    if ((ScalarType(1) - cosom) > std::numeric_limits<ScalarType>::epsilon())
    {
        /* Standard case (slerp) */
        omega = std::acos(cosom);
        sinom = std::sin(omega);
        scale0 = std::sin((ScalarType(1) - t) * omega) / sinom;
        scale1 *= std::sin(t * omega) / sinom;
    }
    else
    {
        /* 'from' and 'to' quaternions are very close, so we can do a linear interpolation */
        scale0 = ScalarType(1) - t;
        scale1 *= t;
    }

    /* Calculate final values */
    return Mix(from, to, scale0, scale1);
}

/**
\brief Returns a smooth 'hermite interpolation' in the range [0, 1].
\remarks This hermite interpolation is: 3x^2 - 2x^3.
*/
template <typename T>
T SmoothStep(const T& x)
{
    return x*x * (T(3) - x*T(2));
}

/**
\brief Returns a smooth 'hermite interpolation' in the range [0, 1].
\remarks This hermite interpolation is: 6x^5 - 15x^4 + 10x^3.
*/
template <typename T>
T SmootherStep(const T& x)
{
    return x*x*x * (x*(x*T(6) - T(15)) + T(10));
}

//! Returns the reciprocal of the specified scalar value.
template <typename T>
T Rcp(const T& x)
{
    return T(1) / x;
}

//! Returns the per-component reciprocal of the specified N-dimensional vector.
template <typename T, std::size_t N>
IVector<T, N> Rcp(const IVector<T, N>& vec)
{
    IVector<T, N> vecRcp;

    for (std::size_t i = 0; i < N; ++i)
        vecRcp[i] = T(1) / vec[i];

    return vecRcp;
}

////! Returns the per-component reciprocal of the specified NxM-dimensional matrix.
//template <typename T, std::size_t N, std::size_t M>
//Matrix<T, N, M> Rcp(const Matrix<T, N, M>& mat)
//{
//    Matrix<T, N, M> matRcp { UninitializeTag{} };

//    for (std::size_t i = 0; i < N*M; ++i)
//        matRcp[i] = T(1) / mat[i];

//    return matRcp;
//}

//! Rescales the specified value 't' from the first range [lower0, upper0] into the second range [lower1, upper1].
template <typename T, typename I>
T Rescale(const T& t, const I& lower0, const I& upper0, const I& lower1, const I& upper1)
{
    /* Return (((t - lower0) / (upper0 - lower0)) * (upper1 - lower1) + lower1) */
    T x = t;
    x -= T(lower0);
    x /= (upper0 - lower0);
    x *= (upper1 - lower1);
    x += T(lower1);
    return x;
}


/* --- Global Operators --- */
}



#endif // IALGEBRA_H
