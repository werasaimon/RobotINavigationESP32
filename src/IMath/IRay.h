/********************************************************************************
 *
 * IRay.h
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

#ifndef IRAY_H
#define IRAY_H

#include "IVector3D.h"

namespace IMath
{

template<class T> class IRay
{
public:

    // -------------------- Attributes -------------------- //

    /// First point of the ray (origin)
    IVector3D<T> Origin;

    /// Second point of the ray
    IVector3D<T> Direction;

    /// Maximum fraction value
    T maxFraction;


    // -------------------- Methods -------------------- //


    /// Constructor default
    SIMD_INLINE IRay(){}

    /// Constructor with arguments
    SIMD_INLINE IRay(const IVector3D<T>& _origin, const IVector3D<T>& _direction, T maxFrac = T(1.0))
        : Origin(_origin),
          Direction(_direction),
          maxFraction(maxFrac)
    {
    }

    /// Copy-constructor
    SIMD_INLINE IRay(const IRay<T>& ray)
        : Origin(ray.Origin),
          Direction(ray.Direction),
          maxFraction(ray.maxFraction)
    {

    }


    /// Overloaded assignment operator
    SIMD_INLINE IRay& operator=(const IRay& ray)
    {
        if (&ray != this)
        {
            Origin = ray.Origin;
            Direction = ray.Direction;
            maxFraction = ray.maxFraction;
        }
        return *this;
    }



    //-------------[ output operator ]------------------------
    /**
        * Output to stream operator
        * @param lhs Left hand side argument of operator (commonly ostream instance).
        * @param rhs Right hand side argument of operator.
        * @return Left hand side argument - the ostream object passed to operator.
        */
    friend std::ostream& operator<<(std::ostream& lhs, const IRay<T>& rhs)
    {
        lhs << "Origin :" << rhs.Origin << " , " <<  "Direction :" << rhs.Direction << " , " << "maxFriction:" << rhs.maxFraction;
        return lhs;
    }

    /**
        * Gets string representation.
        */
    std::string ToString() const
    {
        std::ostringstream oss;
        oss << *this;
        return oss.str();
    }

};




//--------------------------------------
// Typedef shortcuts for IRay
//-------------------------------------

using IRayr    = IRay<Real>;
using IRayf    = IRay<float>;
using IRayd    = IRay<double>;
using IRayi    = IRay<std::int32_t>;
using IRayui   = IRay<std::uint32_t>;
using IRayb    = IRay<std::int8_t>;
using IRayub   = IRay<std::uint8_t>;



} /* namespace */



#endif // IRAY_H
