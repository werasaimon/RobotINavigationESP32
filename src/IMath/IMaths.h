 /********************************************************************************
 *
 * IMaths.h
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

#ifndef MATHS_H
#define MATHS_H

#include "IFunc.h"
#include "IVector2D.h"
#include "IVector3D.h"
#include "IVector4D.h"
#include "ILorentzVector.h"
#include "IMatrix.h"
#include "IMatrix2x2.h"
#include "IMatrix3x3.h"
#include "IMatrix4x4.h"
#include "IQuaternion.h"
#include "IOctonion.h"
#include "IComplex.h"
#include "IRay.h"
#include "ITransform.h"


#include "ISpherical.h"
#include "IPlane.h"


#include "IVector.h"
#include "IAlgebra.h"


#include <limits>

namespace IMath
{

typedef float scalar;

typedef IVector2D<scalar>       Vector2;
typedef IVector3D<scalar>       Vector3;
typedef IVector4D<scalar>       Vector4;
typedef ILorentzVector<scalar>  LorentzVector;
typedef IMatrix2x2<scalar>      Matrix2;
typedef IMatrix3x3<scalar>      Matrix3;
typedef IMatrix4x4<scalar>      Matrix4;
typedef IQuaternion<scalar>     Quaternion;
typedef IRay<scalar>            Ray;
typedef ITransform<scalar>      Transform;
typedef IPlane<scalar>          Plane;

const scalar DECIMAL_SMALLEST = -std::numeric_limits<scalar>::max();
const scalar DECIMAL_LARGEST  =  std::numeric_limits<scalar>::max();



} /* namespace */

#endif // MATHS_H
