 /********************************************************************************
 *
 * ILorentzVector.h
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


#ifndef ILORENTZVECTOR4D_H
#define ILORENTZVECTOR4D_H

#include "IReal.h"
#include "IVector3D.h"


namespace IMath
{



template<class T> class  ILorentzVector
{

    public:


        //! Specifies the typename of the scalar components.
        using ScalarType = T;

        //! Specifies the number of vector components.
        static const std::size_t components = 4;

    //---------------- attribute -------------------//

    // 3 vector component space
    T x;
    T y;
    T z;
    // 1 component time or energy of (x,y,z,t) or (px,py,pz,e)
    T t;


  public:

    // Constructor of the class Vector4D
    SIMD_INLINE ILorentzVector()
    : x(0.0), y(0.0), z(0.0), t(1.0)
    {

    }

    // Constructor with arguments
    SIMD_INLINE ILorentzVector( T newX, T newY, T newZ , T newT )
    : x(newX), y(newY), z(newZ) , t(newT)
    {

    }

    // Copy-constructor
    SIMD_INLINE ILorentzVector(const ILorentzVector<T>& vector)
    : x(vector.x), y(vector.y), z(vector.z) , t(vector.t)
    {

    }


//    // Copy-constructor
//    SIMD_INLINE ILorentzVector(const IVector3D<T>& vector , T time)
//    : x(vector.x), y(vector.y), z(vector.z) , t(time)
//    {

//    }

//    // Copy-constructor
//    SIMD_INLINE ILorentzVector( T time , const IVector3D<T>& vector )
//    : x(vector.x), y(vector.y), z(vector.z) , t(time)
//    {

//    }


    //---------------------- Methods ---------------------//

    SIMD_INLINE void SetToZero()
    {
      x = T(0);
      y = T(0);
      z = T(0);
      t = T(0);
    }


    SIMD_INLINE void SetAllValues(T newX, T newY, T newZ, T newT)
    {
      x = newX;
      y = newY;
      z = newZ;
      t = newT;
    }


    SIMD_INLINE T GetX() const { return x; }
    SIMD_INLINE T GetY() const { return y; }
    SIMD_INLINE T GetZ() const { return z; }
    SIMD_INLINE T GetT() const { return t; }


    SIMD_INLINE T GetP()  const { return (x*x + y*y + z*z); }
    SIMD_INLINE T GetE()  const { return t; }
    SIMD_INLINE T GetEnergy()  const { return t; }


    SIMD_INLINE void SetX(T _x) { x = _x; }
    SIMD_INLINE void SetY(T _y) { y = _y; }
    SIMD_INLINE void SetZ(T _z) { z = _z; }
    SIMD_INLINE void SetW(T _t) { t = _t; }


    SIMD_INLINE void SetXYZM(T x, T y, T z, T m)
    {
        if ( m  >= 0 )
        {
            SetAllValues( x, y, z, ISqrt(x*x+y*y+z*z+m*m) );
        }
        else
        {
            SetAllValues( x, y, z, ISqrt( IMax((x*x+y*y+z*z-m*m), 0. ) ) );
        }
    }


    SIMD_INLINE  void SetPtEtaPhiM(T pt, T eta, T phi, T m)
    {
        pt = Abs(pt);
        SetXYZM(pt*ICos(phi), pt*ISin(phi), pt*ISinh(eta) ,m);
    }

    SIMD_INLINE void SetPtEtaPhiE(T pt, T eta, T phi, T e)
    {
        pt = Abs(pt);
        SetXYZT(pt*ICos(phi), pt*ISin(phi), pt*ISinh(eta) ,e);
    }


    SIMD_INLINE void SetXYZ(const IVector3D<T> &v)
    {
        x = v.x;
        y = v.y;
        z = v.z;
    }

    /// System coordinate . 3 vector component
    SIMD_INLINE IVector3D<T> GetXYZ() const
    {
        return IVector3D<T>(x,y,z);
    }


    ///  Project on space-3D
    SIMD_INLINE IVector3D<T> GetBoostVector() const
    {
       return IVector3D<T>(x/t, y/t, z/t);
    }

    /// Return the square of the length of the vector
    /// Metrices Minkowski Space
    SIMD_INLINE T LengthSquare(const T _SpeedLight = DEFAUL_LIGHT_MAX_VELOCITY_C) const
    {
        const T c = _SpeedLight;
         return (c*c) * (t*t) - (x*x + y*y + z*z);
    }

    /// Return the length of the vector
    SIMD_INLINE T Length(const T _SpeedLight = DEFAUL_LIGHT_MAX_VELOCITY_C) const
    {
        T mm = LengthSquare(_SpeedLight);
        return mm < 0.0 ? -ISqrt(-mm) : ISqrt(mm);
    }


    SIMD_INLINE T GetBeta() const
    {
       return ISqrt(x*x + y*y + z*z) / t;
    }

    SIMD_INLINE T GetGamma() const
    {
       T b = GetBeta();
       return 1.0/ISqrt(1- b*b);
    }



    SIMD_INLINE T Rapidity() const
    {
       //return rapidity
       return 0.5*log( (GetE()+GetZ()) / (GetE()-GetZ()) );
    }



    /**
     * Normalize Unit vector
     */
    SIMD_INLINE ILorentzVector<T> GetUnit() const
    {
        T lengthVector = Length();
        if (lengthVector < MACHINE_EPSILON)
        {
            return *this;
        }
        // Compute and return the unit vector
        T lengthInv = T(1.0) / lengthVector;
        return ILorentzVector<T>( x * lengthInv ,
                                  y * lengthInv ,
                                  z * lengthInv ,
                                  t * lengthInv );
    }


    /**
     * Normalize Unit vector (popular name to methods)
     */
    SIMD_INLINE ILorentzVector<T> Normalized() const
    {
        T lengthVector = Length();
        if (lengthVector < MACHINE_EPSILON)
        {
            return *this;
        }
        // Compute and return the unit vector
        T lengthInv = T(1.0) / lengthVector;
        return ILorentzVector<T>( x * lengthInv ,
                                  y * lengthInv ,
                                  z * lengthInv ,
                                  t * lengthInv );
    }


    /**
    * Inverse vector
    */
    SIMD_INLINE ILorentzVector<T> GetInverse() const
    {
        return ILorentzVector<T>( T(1.0/x) , T(1.0/y) , T(1.0/z) , T(1.0/t));
    }





    //---------------[ vector aritmetic operator ]--------------
    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T> operator+(const ILorentzVector<T>& rhs) const
    {
        return ILorentzVector<T>(x + rhs.x, y + rhs.y, z + rhs.z, t + rhs.t);
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T> operator-(const ILorentzVector<T>& rhs) const
    {
        return ILorentzVector<T>(x - rhs.x, y - rhs.y, z - rhs.z, t - rhs.t);
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T> operator*(const ILorentzVector<T> rhs) const
    {
        return ILorentzVector<T>(x * rhs.x, y * rhs.y, z * rhs.z, t * rhs.t);
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T> operator/(const ILorentzVector<T>& rhs) const
    {
        return ILorentzVector<T>(x / rhs.x, y / rhs.y, z / rhs.z, t / rhs.t);
    }

    //------------------------------ Friends ----------------------------------------//

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    friend SIMD_INLINE ILorentzVector<T> operator*(T number, const ILorentzVector<T>& vector)
    {
        return ILorentzVector<T>(number * vector.x, number * vector.y, number * vector.z , number * vector.t);
    }


    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    friend SIMD_INLINE  ILorentzVector<T> operator/( T number , const ILorentzVector<T>& vector )
    {
        return ILorentzVector<T>(vector.x / number, vector.y / number, vector.z / number , vector.t / number);
    }

    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T>& operator+=(const ILorentzVector<T>& rhs)
    {
    	x += rhs.x;
    	y += rhs.y;
    	z += rhs.z;
    	t += rhs.t;
    	return *this;
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T>& operator-=(const ILorentzVector<T>& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        t -= rhs.t;
        return *this;
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T>& operator*=(const ILorentzVector<T>& rhs)
    {
        x *= rhs.x;
        y *= rhs.y;
        z *= rhs.z;
        t *= rhs.t;
        return *this;
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T>& operator/=(const ILorentzVector<T>& rhs)
    {
        x /= rhs.x;
        y /= rhs.y;
        z /= rhs.z;
        t /= rhs.t;
        return *this;
    }


    //--------------[ equiality operator ]------------------------
    /**
     * Equality test operator
     * @param rhs Right hand side argument of binary operator.
     * @note Test of equality is based of threshold MACHINE_EPSILON value. To be two
     * values equal, must satisfy this condition | lhs.x - rhs.y | < MACHINE_EPSILON,
     * same for y-coordinate, z-coordinate, and w-coordinate.
     */
    SIMD_INLINE bool operator==(const ILorentzVector<T>& rhs) const
    {
        return IAbs(x - rhs.x) < MACHINE_EPSILON &&
               IAbs(y - rhs.y) < MACHINE_EPSILON &&
               IAbs(z - rhs.z) < MACHINE_EPSILON &&
               IAbs(t - rhs.t) < MACHINE_EPSILON;
    }

    /**
     * Inequality test operator
     * @param rhs Right hand side argument of binary operator.
     * @return not (lhs == rhs) :-P
     */
    SIMD_INLINE bool operator!=(const ILorentzVector<T>& rhs) const
    {
    	return !(*this == rhs);
    }


    //----------------[ access operators ]-------------------
     /**
      * Copy operator
      * @param rhs Right hand side argument of binary operator.
      */
    SIMD_INLINE ILorentzVector<T> operator=(const ILorentzVector<T>& rhs)
     {
         x = rhs.x;
         y = rhs.y;
         z = rhs.z;
         t = rhs.t;
         return *this;
     }

     /**
      * Copy casting operator
      * @param rhs Right hand side argument of binary operator.
      */
     template<class FromT>
     SIMD_INLINE ILorentzVector<T> operator=(const ILorentzVector<FromT>& rhs)
     {
         x = static_cast<T>(rhs.x);
         y = static_cast<T>(rhs.y);
         z = static_cast<T>(rhs.z);
         t = static_cast<T>(rhs.t);
         return *this;
     }


    /**
     * Array access operator
     * @param n Array index
     * @return For n = 0, reference to x coordinate, n = 1
     * reference to y coordinate, n = 2 reference to z,
     * else reference to t-time coordinate.
     */
     SIMD_INLINE T & operator[](int n)
    {
        static_assert(sizeof(*this) == sizeof(T[components]), "");
    	assert(n >= 0 && n < 4);
    	return (&x)[n];
    }

    /**
     * Array access operator
     * @param n Array index
     * @return For n = 0, reference to x coordinate, n = 1
     * reference to y coordinate, n = 2 reference to z,
     * else reference to t time-coordinate.
     */
    SIMD_INLINE const T & operator[](int n) const
    {
        static_assert(sizeof(*this) == sizeof(T[components]), "");
    	assert(n >= 0 && n < 4);
    	return (&x)[n];
    }

    //-------------[ unary operations ]--------------------------
    /**
     * Unary negate operator
     * @return negated vector
     */
    SIMD_INLINE ILorentzVector<T> operator-() const
    {
        return ILorentzVector<T>(-x, -y, -z, -t);
    }

    //--------------[ scalar vector operator ]--------------------

    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T> operator+(T rhs) const
    {
        return ILorentzVector<T>(x + rhs, y + rhs, z + rhs, t + rhs);
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T> operator-(T rhs) const
    {
        return ILorentzVector<T>(x - rhs, y - rhs, z - rhs, t - rhs);
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T> operator*(T rhs) const
    {
        return ILorentzVector<T>(x * rhs, y * rhs, z * rhs, t * rhs);
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T> operator/(T rhs) const
    {
        return ILorentzVector<T>(x / rhs, y / rhs, z / rhs, t / rhs);
    }

    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T>& operator+=(T rhs)
    {
        x += rhs;
        y += rhs;
        z += rhs;
        t += rhs;
        return *this;
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T>& operator-=(T rhs)
    {
        x -= rhs;
        y -= rhs;
        z -= rhs;
        t -= rhs;
        return *this;
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T>& operator*=(T rhs)
    {
        x *= rhs;
        y *= rhs;
        z *= rhs;
        t *= rhs;
        return *this;
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T>& operator/=(T rhs)
    {
        x /= rhs;
        y /= rhs;
        z /= rhs;
        t /= rhs;
        return *this;
    }


    /**
     * Dot product of two vectors.
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE T Dot(const ILorentzVector<T>& rhs) const
    {
         return t*rhs.t - z*rhs.z - y*rhs.y - x*rhs.x;
    }

    /**
     * Cross tri product operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE ILorentzVector<T> Cross(const ILorentzVector<T>& b , const ILorentzVector<T>& c) const
    {

        //Precompute some 2x2 matrix determinants for speed
         T Pxy = b.x*c.y - c.x*b.y;
         T Pxz = b.x*c.z - c.x*b.z;
         T Pxw = b.x*c.t - c.x*b.t;
         T Pyz = b.y*c.z - c.y*b.z;
         T Pyw = b.y*c.t - c.y*b.t;
         T Pzw = b.z*c.t - c.z*b.t;

          return ILorentzVector<T>
          (
             y*Pzw - z*Pyw + t*Pyz,    //Note the lack of 'x' in this line
             z*Pxw - x*Pzw - t*Pxz,    //y, Etc.
             x*Pyw - y*Pxw + t*Pxy,
             y*Pxz - x*Pyz - z*Pxy
          );
    }




    //================================  Method Gerglocema =======================================//

//    SIMD_INLINE void Boost(T bx, T by, T bz , const T _SpeedLight = DEFAUL_LIGHT_MAX_VELOCITY_C)
//    {
//        /**
//       //Boost this Lorentz vector
//       T b2 = bx*bx + by*by + bz*bz;
//       T gamma = 1.0 / ISqrt(1.0 - b2);
//       T bp = bx*x + by*y + bz*z;
//       T gamma2 = b2 > 0 ? (gamma - 1.0)/b2 : 0.0;

//       SetX(GetX() + gamma2*bp*bx + gamma*bx*GetT());
//       SetY(GetY() + gamma2*bp*by + gamma*by*GetT());
//       SetZ(GetZ() + gamma2*bp*bz + gamma*bz*GetT());
//       SetT(gamma*(GetT() + bp));
//       **/

//        /// method Gerglocema
//        *this = CreateGerglocemaBoost( *this , IVector3D<T>(bx,by,bz) , _SpeedLight);
//    }

//    SIMD_INLINE void Boost(T bx, T by, T bz , T gamma , const T _SpeedLight = DEFAUL_LIGHT_MAX_VELOCITY_C)
//    {
//    	/// method Gerglocema
//        *this = CreateGerglocemaBoost( *this , IVector3D<T>(bx,by,bz) , gamma , _SpeedLight);
//    }

    SIMD_INLINE void GerglocemaBoostInvert( const IVector3D<T> &b , const T _SpeedLight = DEFAUL_LIGHT_MAX_VELOCITY_C)
    {
        /// Invert method Gerglocema
        *this = CreateGerglocemaBoostInvert( *this , b , _SpeedLight );
    }

    SIMD_INLINE void GerglocemaBoost( T gamma , const IVector3D<T> &b ,  const T _SpeedLight = DEFAUL_LIGHT_MAX_VELOCITY_C)
    {
        /// Method Gerglocema
        *this = CreateGerglocemaBoost( gamma  , *this , b , _SpeedLight );
    }


    //===========================================================================================//




     //================================ Plugin =====================================//

    /**
     *  Method Gerglocema Invert-Method ^ -1
     */
    static SIMD_INLINE ILorentzVector<T> CreateGerglocemaBoostInvert( const ILorentzVector<T> &pos , const IVector3D<T> &v , const T _SpeedLight = DEFAUL_LIGHT_MAX_VELOCITY_C)
    {
        /// Light speed
        const T c = _SpeedLight;

        /// Invert gamma factor
        T gamma = 1.0 * ISqrt( 1.0 + (v.Dot(v)) / (c*c) );

        /// method Gerglocema
        ILorentzVector<T> res;
        IVector3D<T> ortogonalPredikat = (pos.GetXYZ().Cross(v)).Cross(v);
        res.SetXYZ(((pos.GetXYZ() + v * pos.t ) / gamma) - (T(1)/(c*c))*(T(1)/(gamma*(1+gamma))) * ortogonalPredikat);
        res.t =  (pos.t + (v.Dot(pos.GetXYZ())/(c*c))) / gamma ;
        return res;
    }


    /**
     *  Method Gerglocema
     */
    static SIMD_INLINE ILorentzVector<T> CreateGerglocemaBoost( T gamma , const ILorentzVector<T> &pos , const IVector3D<T> &v , const T _SpeedLight = DEFAUL_LIGHT_MAX_VELOCITY_C)
    {
        ///Light speed
        const T c = _SpeedLight;

        /// method Gerglocema
        ILorentzVector<T> res;
        IVector3D<T> ortogonalPredikat = (pos.GetXYZ().Cross(v)).Cross(v);

        res.x = ((pos.x + v.x * pos.t ) / gamma) + (T(1)/(c*c))*(T(1)/(gamma*(1+gamma))) * ortogonalPredikat.x;
        res.y = ((pos.y + v.y * pos.t ) / gamma) + (T(1)/(c*c))*(T(1)/(gamma*(1+gamma))) * ortogonalPredikat.y;
        res.z = ((pos.z + v.z * pos.t ) / gamma) + (T(1)/(c*c))*(T(1)/(gamma*(1+gamma))) * ortogonalPredikat.z;
        res.t =  (pos.t + (v.Dot(pos.GetXYZ())/(c*c))) / gamma ;

        return res;
    }


    //-------------[ output operator ]------------------------
     /**
     * Output to stream operator
     * @param lhs Left hand side argument of operator (commonly ostream instance).
     * @param rhs Right hand side argument of operator.
     * @return Left hand side argument - the ostream object passed to operator.
     */
     friend std::ostream& operator<<(std::ostream& lhs, const ILorentzVector<T>& rhs)
     {
         lhs << "[" << rhs[0] << "," << rhs[1] << "," << rhs[2] << "," << rhs[3] << "]";
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


template<class T> const
static ILorentzVector<T> Cross(const ILorentzVector<T>& a, const ILorentzVector<T>& b , const ILorentzVector<T>& c)
{
    return a.Cross(b,c);
}

template<class T> const
static T Dot(const ILorentzVector<T>& a, const ILorentzVector<T>& b)
{
    return a.Dot(b);
}


//--------------------------------------
// Typedef shortcuts for 4D LorentzVector
//-------------------------------------

using ILorentzVectorr    = ILorentzVector<Real>;
using ILorentzVectorf    = ILorentzVector<float>;
using ILorentzVectord    = ILorentzVector<double>;
using ILorentzVectori    = ILorentzVector<std::int32_t>;
using ILorentzVectorui   = ILorentzVector<std::uint32_t>;
using ILorentzVectorb    = ILorentzVector<std::int8_t>;
using ILorentzVectorub   = ILorentzVector<std::uint8_t>;



} /* namespace */

#endif // ILORENTZVECTOR4D_H
