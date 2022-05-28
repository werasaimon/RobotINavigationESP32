 /********************************************************************************
 *
 * IFunc.h
 *
 * IMath : 3d_math library,
 * Copyright (c)  *
 * Created on: 3 July. 2018 г.
 * Author: werasaimon                                     *
 *********************************************************************************
 *                                                                               *
 * This software is provided 'as-is', without any express or implied warranty.   *
 * In no event will the authors be held liable for any damages arising from the  *
 * use of this software.                                                         *
 *                                                                               *
 * Permission is granted to anyone to use this software for any purpose,         *
 * including commercial applications, and to alter it and redistribute it        *
 * freely, subject to the following restrictions:                                *
 *                                                                               *
 * 1. The origin of this software must not be misrepresented; you must not claim *
 *    that you wrote the original software. If you use this software in a        *
 *    product, an acknowledgment in the product documentation would be           *
 *    appreciated but is not required.                                           *
 *                                                                               *
 * 2. Altered source versions must be plainly marked as such, and must not be    *
 *    misrepresented as being the original software.                             *
 *                                                                               *
 * 3. This notice may not be removed or altered from any source distribution.    *
 *                                                                               *
 ********************************************************************************/


#ifndef FUNC_H
#define FUNC_H


#ifdef _DEBUG
#include <crtdbg.h>
#define _CRTDBG_MAP_ALLOC
#endif

//Liberies
//#include <cmath>

#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <cassert>
#include <cstdlib>
#include <limits>
#include <math.h>
#include <cmath>


namespace IMath
{

#ifndef M_PI
#define M_PI  acos(-1) // 3.14159265358979323846  /* pi */
#endif



/// Light Velocity c = 300000.kilometers / 1.second
#define DEFAUL_LIGHT_MAX_VELOCITY_C 1.0

//-------------------------------------------------------------------------------
//-- Typedefs, Structs ----------------------------------------------------------
//-------------------------------------------------------------------------------


#ifndef SIMD_INLINE
# if defined(_MSC_VER) && (_MSC_VER >= 1200)
#  define SIMD_INLINE __forceinline
# else
#  define SIMD_INLINE __inline
# endif
#endif


#define ATTRIBUTE_ALIGNED16(a)  __declspec(align(16)) a
#define ATTRIBUTE_ALIGNED64(a)  __declspec(align(64)) a
#define ATTRIBUTE_ALIGNED128(a) __declspec(align(128)) a



#define FLT_MIN 1.175494351e-38F /* min positive value */
#define FLT_MAX 3.402823466e+38F /* max value */


#define FLT_EPSILON 1.1920928955078125E-7f
#define MACHINE_EPSILON 10E-7f


// FLOATING POINT VALIDITY
template<typename T> bool is_infinite(const T &value)
{
    T max_value = (std::numeric_limits<T>::max)();
    T min_value = -max_value;
    return !(min_value <= value && value <= max_value);
}

template<typename T> bool is_nan(const T &value)
{
    // True if NAN
    return value != value;
}

template<typename T> bool is_valid(const T &value)
{
    return !is_infinite(value) && !is_nan(value);
}


    //-----------------------------------------------------//

    /**
    \brief Computes a linear interpolation between the point 'a' and the point 'b'.
    \return Equivalent to: a*(1-t) + b*t
    */
    template <typename T, typename I>
    void Lerp(T& x, const T& a, const T& b, const I& t)
    {
        x = b;
        x -= a;
        x *= t;
        x += a;
    }

    /**
    \brief Computes a linear interpolation between the point 'a' and the point 'b'.
    \return Equivalent to: a*(1-t) + b*t
    */
    template <typename T, typename I>
    T Lerp(const T& a, const T& b, const I& t)
    {
        /* Return (b - a) * t + a */
        T x = b;
        x -= a;
        x *= t;
        x += a;
        return x;
    }

    //-----------------------------------------------------//



// ---------- Mathematics functions ---------- //

template<typename T> SIMD_INLINE T IMin(T a, T b)
{
    return (a > b) ? b : a;
}

template<typename T> SIMD_INLINE T IMax(T a, T b)
{
    return (a < b) ? b : a;
}

template<typename T> SIMD_INLINE T IMin(T a, T b, T c)
{
    return IMin<T>(IMin<T>(a, b), c);
}

template<typename T> SIMD_INLINE T IMax(T a, T b, T c)
{
    return IMax<T>(IMax<T>(a, b), c);
}

template<typename T> SIMD_INLINE T IClamp(T a, T min, T max)
{
    return IMax<T>(IMin<T>(a, max), min);
}

template<typename T> SIMD_INLINE T IWrap(T a, T min, T max)
{
    return (a < min) ? max - (min - a) : (a > max) ? min - (max - a) : a;
}

template<typename T> SIMD_INLINE void ISwap(T& a, T& b)
{
    T c = a;
    a = b;
    b = c;
}

//-----------------------------------------------------//


template<typename T> SIMD_INLINE T ISign(T  x)
{
    return (x < 0.0f) ? -1.0f : 1.0f;
}


template<typename T> SIMD_INLINE T IPi(void)
{
    static const T  gPi = atan(1.0f) * 4.0f;
    return gPi;
}
template<typename T> SIMD_INLINE T ITwoPi(void)
{
    static const T  gTwoPi = atan(1.0f) * 8.0f;
    return gTwoPi;
}
template<typename T> SIMD_INLINE T IModulo(T  x, T  div)
{
    return  fmod( x, div);
}
template<typename T> SIMD_INLINE T IAbs(T  x)
{
    return fabs(x);
}

template<typename T>  bool IsZero( T a, T epsilon = MACHINE_EPSILON )
{
    return (IAbs(a) <= epsilon);
}
//-------------------------------------------------------------------------------
// Is this floating point value close to zero?
//-------------------------------------------------------------------------------
template<typename T>  bool IIsZero( T a, T epsilon = MACHINE_EPSILON )
{
    return (IAbs(a) <= epsilon);
}


template<typename T> SIMD_INLINE T IDegreesToRadians(T  Degrees)
{
    return Degrees * (M_PI / 180.0f);
}
template<typename T> SIMD_INLINE T IRadiansToDegrees(T  Radians)
{
    return Radians * (180.0f / M_PI);
}


template<typename T> SIMD_INLINE T ISin(T  Radians)
{
    return  sin(Radians);
}
template<typename T> SIMD_INLINE T ICos(T  Radians)
{
    return  cos(Radians);
}

template<typename T> SIMD_INLINE T ISinh(T  Radians)
{
    return  sinh(Radians);
}
template<typename T> SIMD_INLINE T ICosh(T  Radians)
{
    return  cosh(Radians);
}

template<typename T> SIMD_INLINE T ISinc_pi(T x)
{
  return  sin(x) / x;
}

template<typename T> SIMD_INLINE T ICosc_pi(T x)
{
  return  cos(x) / x;
}

template<typename T> SIMD_INLINE T ISinhc_pi(T x)
{
    return  sinh(x) / x;
}

template<typename T> SIMD_INLINE T ICoshc_pi(T x)
{
    return  cosh(x) / x;
}


template<typename T> SIMD_INLINE T ITan(T  Radians)
{
    return  tan(Radians);
}
template<typename T> SIMD_INLINE T IAtan(T  X)
{
    return  atan(X);
}

template<typename T> SIMD_INLINE T IASin(T  X)
{
    return  asin(X);
}
template<typename T> SIMD_INLINE T IACos(T  X)
{
    return  acos(X);
}


template<typename T> SIMD_INLINE T IAtan2(T  X , T  Y)
{
    return  atan2(X , Y);
}

template<typename T> SIMD_INLINE T ILog(T  X)
{
    return  log(X);
}


template<typename T> SIMD_INLINE T ILog2(T  X )
{
    return  log2(X);
}

template<typename T> SIMD_INLINE T IExp(T  X )
{
    return  exp(X);
}

template<typename T> SIMD_INLINE T IExp2(T  X )
{
    return  exp2(X);
}

template<typename T> SIMD_INLINE T ISqrt(T  In)
{
    return sqrt(In);
}

template<typename T> SIMD_INLINE T IRand(T  r = 1.0f)
{
    return rand() / ( RAND_MAX) * r;
}
template<typename T> SIMD_INLINE T IRand(T  min, T  max)
{
    return min + Rand(max - min);
}

template<typename T> SIMD_INLINE T IPow(T  X , T  Y)
{
    return  pow(X , Y);
}



template<typename T> SIMD_INLINE T IMin3(T a, T b, T c)
{
    return IMin(IMin(a, b), c);
}


template<typename T> SIMD_INLINE T IMax3(T a, T b, T c)
{
    return IMax(IMax(a, b), c);
}


/// Function to test if two real numbers are (almost) equal
/// We test if two numbers a and b are such that (a-b) are in [-EPSILON; EPSILON]
template<typename T> SIMD_INLINE bool IApproxEqual(T a, T b, T epsilon = MACHINE_EPSILON)
{
    return (IAbs(a - b) < T(epsilon) );
}


template<typename T> SIMD_INLINE bool ISameSign(T a, T b)
{
    return a * b >= T(0.0);
}


}  /* namespace */


#endif // FUNC_H
