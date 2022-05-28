/********************************************************************************
 *
 * IVector4D.h
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

#ifndef RPVECTOR4D_H
#define RPVECTOR4D_H

#include "IReal.h"
#include "IVector3D.h"


namespace IMath
{
/**
 * Class for four dimensional vector.
  * There are four ways of accessing vector components.
 * Let's have <code>Vector4f v</code>, you can either:
 * <ul>
 * 	<li>access as position in projective space (x,y,z,w) &mdash; <code>v.x = v.y = v.z = v.w = 1;</code></li>
 * 	<li>access as texture coordinate (s,t,u,v) &mdash; <code>v.s = v.t = v.u = v.v = 1;</code></li>
 * 	<li>access as color (r,g,b,a) &mdash; <code>v.r = v.g = v.b = v.a = 1;</code></li>
 * 	<li>access via operator[] &mdash; <code>v[0] = v[1] = v[2] = v[3] = 1;</code></li>
 * </ul>
 */
template<class T>
class IVector4D
{

public:


    //! Specifies the typename of the scalar components.
    using ScalarType = T;

    //! Specifies the number of vector components.
    static const std::size_t components = 4;

    //-------------------- Attributes --------------------//

    union
    {
        /**
         * First element of vector, alias for R-coordinate.
         * For color notation.
         */
        T r
        /**
         * First element of vector, alias for X-coordinate.
         */;
        T x;


    };

    union
    {
        /**
         * Second element of vector, alias for G-coordinate.
         * For color notation.
         */
        T g;
        /**
         * Second element of vector, alias for Y-coordinate.
         */
        T y;

    };

    union
    {
        /**
         * Third element of vector, alias for B-coordinate.
         * For color notation.
         */
        T b;
        /**
         * Third element of vector, alias for Z-coordinate.
         */
        T z;

    };

    union
    {
        /**
         * First element of vector, alias for W-coordinate.
         * @note For vectors (such as normals) should be Set to 0.0
         * For vertices should be Set to 1.0
         */
        T w;
        /**
         * Fourth element of vector, alias for A-coordinate.
         * For color notation. This represnt aplha chanell
         */
        T a;
    };

public:

    //----------------[ constructors ]--------------------------
    /**
     * Creates and Sets to (0,0,0,0)
     */
    SIMD_INLINE IVector4D()
        : x(0), y(0), z(0), w(0)
    {
    }




    SIMD_INLINE explicit IVector4D(const T& scalar)
        : x( scalar ) ,
          y( scalar ) ,
          z( scalar ) ,
          w( scalar )
    {

    }

    /**
     * Creates and Sets to (x,y,z,w)
     * @param nx initial x-coordinate value (R)
     * @param ny initial y-coordinate value (G)
     * @param nz initial z-coordinate value (B)
     * @param nw initial w-coordinate value (Alpha)
     */
    SIMD_INLINE IVector4D(T nx, T ny, T nz, T nw)
        : x(nx), y(ny), z(nz), w(nw)
    {
    }


    /**
     * Creates and Sets to (n,w)
     * @param nx initial n.x-coordinate value (R)
     * @param ny initial n.y-coordinate value (G)
     * @param nz initial n.z-coordinate value (B)
     * @param nw initial w-coordinate value (Alpha)
     */
    SIMD_INLINE IVector4D(const IVector3D<T>& n , T nw)
        : x(n.x), y(n.y), z(n.z), w(nw)
    {
    }

    /**
     * Copy constructor.
     * @param src Source of data for new created IVector4D instance.
     */
    SIMD_INLINE IVector4D(const IVector4D<T>& src)
        : x(src.x), y(src.y), z(src.z), w(src.w)
    {
    }

    /**
     * Copy casting constructor.
     * @param src Source of data for new created IVector4D instance.
     */
    template<class FromT>
    SIMD_INLINE IVector4D(const IVector4D<FromT>& src)
        : x(static_cast<T>(src.x)),
          y(static_cast<T>(src.y)),
          z(static_cast<T>(src.z)),
          w(static_cast<T>(src.w))
    {
    }


    //---------------------- Methods ---------------------//

    SIMD_INLINE void SetToZero()
    {
        x = T(0);
        y = T(0);
        z = T(0);
        w = T(0);
    }


    SIMD_INLINE void SetAllValues(T newX, T newY, T newZ, T newW)
    {
        x = newX;
        y = newY;
        z = newZ;
        w = newW;
    }


    SIMD_INLINE T GetX() const { return x; }
    SIMD_INLINE T GetY() const { return y; }
    SIMD_INLINE T GetZ() const { return z; }
    SIMD_INLINE T GetW() const { return w; }

    SIMD_INLINE IVector2D<T> GetXY()  const { return IVector2D<T>(x,y); }
    SIMD_INLINE IVector3D<T> GetXYZ() const { return IVector3D<T>(x,y,z); }


    SIMD_INLINE void SetX(T _x) { x = _x; }
    SIMD_INLINE void SetY(T _y) { y = _y; }
    SIMD_INLINE void SetZ(T _z) { z = _z; }
    SIMD_INLINE void SetW(T _w) { w = _w; }


    SIMD_INLINE void SetXY(const IVector2D<T>& _v)
    {
        x=_v.x;
        y=_v.y;
    }

    SIMD_INLINE void SetXYZ(const IVector3D<T>& _v)
    {
        x=_v.x;
        y=_v.y;
        z=_v.z;
    }


    //----------------[ access operators ]-------------------
    /**
     * Copy operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector4D<T> operator=(const IVector4D<T>& rhs)
    {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        w = rhs.w;
        return *this;
    }

    /**
     * Copy casting operator
     * @param rhs Right hand side argument of binary operator.
     */
    template<class FromT>
    SIMD_INLINE IVector4D<T> operator=(const IVector4D<FromT>& rhs)
    {
        x = static_cast<T>(rhs.x);
        y = static_cast<T>(rhs.y);
        z = static_cast<T>(rhs.z);
        w = static_cast<T>(rhs.w);
        return *this;
    }

    /**
     * Array access operator
     * @param n Array index
     * @return For n = 0, reference to x coordinate, n = 1
     * reference to y coordinate, n = 2 reference to z,
     * else reference to w coordinate.
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
     * else reference to w coordinate.
     */
    SIMD_INLINE const T & operator[](int n) const
    {
        static_assert(sizeof(*this) == sizeof(T[components]), "");
        assert(n >= 0 && n < 4);
        return (&x)[n];
    }

    //---------------[ vector aritmetic operator ]--------------
    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector4D<T> operator+(const IVector4D<T>& rhs) const
    {
        return IVector4D<T>(x + rhs.x, y + rhs.y, z + rhs.z, w + rhs.w);
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector4D<T> operator-(const IVector4D<T>& rhs) const
    {
        return IVector4D<T>(x - rhs.x, y - rhs.y, z - rhs.z, w - rhs.w);
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector4D<T> operator*(const IVector4D<T> rhs) const
    {
        return IVector4D<T>(x * rhs.x, y * rhs.y, z * rhs.z, w * rhs.w);
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector4D<T> operator/(const IVector4D<T>& rhs) const
    {
        return IVector4D<T>(x / rhs.x, y / rhs.y, z / rhs.z, w / rhs.w);
    }

    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector4D<T>& operator+=(const IVector4D<T>& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        w += rhs.w;
        return *this;
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector4D<T>& operator-=(const IVector4D<T>& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        w -= rhs.w;
        return *this;
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector4D<T>& operator*=(const IVector4D<T>& rhs)
    {
        x *= rhs.x;
        y *= rhs.y;
        z *= rhs.z;
        w *= rhs.w;
        return *this;
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector4D<T>& operator/=(const IVector4D<T>& rhs)
    {
        x /= rhs.x;
        y /= rhs.y;
        z /= rhs.z;
        w /= rhs.w;
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
    SIMD_INLINE bool operator==(const IVector4D<T>& rhs) const
    {
        return IAbs(x - rhs.x) < MACHINE_EPSILON &&
               IAbs(y - rhs.y) < MACHINE_EPSILON &&
               IAbs(z - rhs.z) < MACHINE_EPSILON &&
               IAbs(w - rhs.w) < MACHINE_EPSILON;
    }

    /**
     * Inequality test operator
     * @param rhs Right hand side argument of binary operator.
     * @return not (lhs == rhs) :-P
     */
    SIMD_INLINE bool operator!=(const IVector4D<T>& rhs) const
    {
        return !(*this == rhs);
    }

    //-------------[ unary operations ]--------------------------
    /**
     * Unary negate operator
     * @return negated vector
     */
    SIMD_INLINE IVector4D<T> operator-() const
    {
        return IVector4D<T>(-x, -y, -z, -w);
    }

    //--------------[ scalar vector operator ]--------------------

    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector4D<T> operator+(T rhs) const
    {
        return IVector4D<T>(x + rhs, y + rhs, z + rhs, w + rhs);
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector4D<T> operator-(T rhs) const
    {
        return IVector4D<T>(x - rhs, y - rhs, z - rhs, w - rhs);
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector4D<T> operator*(T rhs) const
    {
        return IVector4D<T>(x * rhs, y * rhs, z * rhs, w * rhs);
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector4D<T> operator/(T rhs) const
    {
        return IVector4D<T>(x / rhs, y / rhs, z / rhs, w / rhs);
    }

    //------------------------------ Friends ----------------------------------------//

    /**
      * Multiplication operator
      * @param rhs Right hand side argument of binary operator.
      */
    friend SIMD_INLINE IVector4D<T> operator*(T number, const IVector4D<T>& vector)
    {
        return IVector4D<T>(number * vector.x, number * vector.y, number * vector.z , number * vector.w);
    }


    /**
      * Division operator
      * @param rhs Right hand side argument of binary operator.
      */
    friend SIMD_INLINE  IVector4D<T> operator/( T number , const IVector4D<T>& vector )
    {
        return IVector4D<T>(vector.x / number, vector.y / number, vector.z / number , vector.w / number);
    }


    //------------------------------ Dynamics ----------------------------------------//

    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector4D<T>& operator+=(T rhs)
    {
        x += rhs;
        y += rhs;
        z += rhs;
        w += rhs;
        return *this;
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector4D<T>& operator-=(T rhs)
    {
        x -= rhs;
        y -= rhs;
        z -= rhs;
        w -= rhs;
        return *this;
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector4D<T>& operator*=(T rhs)
    {
        x *= rhs;
        y *= rhs;
        z *= rhs;
        w *= rhs;
        return *this;
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector4D<T>& operator/=(T rhs)
    {
        x /= rhs;
        y /= rhs;
        z /= rhs;
        w /= rhs;
        return *this;
    }

    //------------------------------ Dot . Cross ----------------------------------------//

    /**
     * Dot product of two vectors.
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE T Dot(const IVector4D<T>& rhs) const
    {
        return x * rhs.x + y * rhs.y + z * rhs.z + w * rhs.w;
    }

    /**
     * Cross tri product operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector4D<T> Cross(const IVector4D<T>& b , const IVector4D<T>& c) const
    {

        //Precompute some 2x2 matrix determinants for speed
        T Pxy = b.x*c.y - c.x*b.y;
        T Pxz = b.x*c.z - c.x*b.z;
        T Pxw = b.x*c.w - c.x*b.w;
        T Pyz = b.y*c.z - c.y*b.z;
        T Pyw = b.y*c.w - c.y*b.w;
        T Pzw = b.z*c.w - c.z*b.w;

        return IVector4D<T>
                (
                    y*Pzw - z*Pyw + w*Pyz,    //Note the lack of 'x' in this line
                    z*Pxw - x*Pzw - w*Pxz,    //y, Etc.
                    x*Pyw - y*Pxw + w*Pxy,
                    y*Pxz - x*Pyz - z*Pxy
                    );
    }


    //-------------[ size operations ]---------------------------


    /**
     * Return square of length.
     * @return length ^ 2
     * @note This method is faster then length(). For comparison
     * of length of two vector can be used just this value, instead
     * of more expensive length() method.
     */
    SIMD_INLINE T LengthSquare() const
    {
        return x * x + y * y + z * z + w * w;
    }


    /**
     * Get length of vector.
     * @return lenght of vector
     */
    SIMD_INLINE T Length() const
    {
        return ISqrt(x * x + y * y + z * z + w * w);
    }

    /**
     * Normalize vector
     */
    SIMD_INLINE void Normalize()
    {
        T s = Length();
        x /= s;
        y /= s;
        z /= s;
        w /= s;
    }

    /**
     * Normalize Unit vector
     */
    SIMD_INLINE IVector4D<T> GetUnit() const
    {
        T lengthVector = Length();
        if (lengthVector < MACHINE_EPSILON)
        {
            return *this;
        }
        // Compute and return the unit vector
        T lengthInv = T(1.0) / lengthVector;
        return IVector4D<T>( x * lengthInv,
                             y * lengthInv,
                             z * lengthInv,
                             w * lengthInv);
    }

    /**
     * Normalize Unit vector (popular name to methods)
     */
    SIMD_INLINE IVector4D<T>  Normalized() const
    {
        T lengthVector = Length();
        if (lengthVector < MACHINE_EPSILON)
        {
            return *this;
        }
        // Compute and return the unit vector
        T lengthInv = T(1.0) / lengthVector;
        return IVector4D<T>( x * lengthInv,
                             y * lengthInv,
                             z * lengthInv,
                             w * lengthInv);
    }

    /**
    * Inverse vector
    */
    SIMD_INLINE IVector4D<T> Inverse() const
    {
        return IVector4D<T>( T(1.0/x) , T(1.0/y) , T(1.0/z) , T(1.0/w));
    }


    /**
    * Transpose Vector
    */
    SIMD_INLINE IVector4D<T> Transpose() const
    {
          return IVector4D<T>(w,z,y,x);
    }



    /**
    * @Randomize Vector
    */
    SIMD_INLINE T Random( T l, T h ) const
    {
      T a = (T)rand( );
      a /= RAND_MAX;
      a = (h - l) * a + l;
      return a;
    }



    //    //--------------[ Logic-Plane ]-----------------------


    //    /**
    //    * Build make plane
    //    */
    //    SIMD_INLINE IVector4D<T> BuildPlan(const IVector3D<T> & p_point1, const IVector3D<T> &p_normal)
    //    {
    //       IVector3D<T> normal;
    //       normal = (p_normal.normalized());
    //       w  = normal.dot(p_point1);
    //       x  = normal.x;
    //       y  = normal.y;
    //       z  = normal.z;
    //       return *this;
    //    }


    //    /**
    //    * Intersection plane to direction
    //    */
    //    SIMD_INLINE bool RayInter( IVector3D<T> &interPoint, const IVector3D<T> &position, const IVector3D<T> &direction)
    //    {
    //        float den = IVector3D<T>(x,y,z).dot(direction);

    //        IVector3D<T> P = position;

    ///*        if( IAbs(den) < 0.00001 )
    //        {
    //            P = IVector3D<T>(x,y,z);
    //        }
    //        else*/ if( IAbs(den) < T(0.000001) )
    //        {
    //           return false;
    //        }


    //        IVector3D<T> tmp = (IVector3D<T>(x,y,z) * w) - P;
    //        interPoint = P + (IVector3D<T>(x,y,z).dot(tmp) / den) * direction;

    //        return true;
    //    }


    SIMD_INLINE IVector3D<T> ClosestPoint( const IVector3D<T>& point ) const
    {
        IVector3D<T> Normal(x,y,z);
        return point - (Normal.Dot(point) - w)*Normal;
    }


    //--------------[ misc. operations ]-----------------------

    /**
     * Linear interpolation of two vectors
     * @param fact Factor of interpolation. For translation from position
     * of this vector to vector r, values of factor goes from 0.0 to 1.0.
     * @param r Second Vector for interpolation
     * @note However values of fact parameter are reasonable only in interval
     * [0.0 , 1.0], you can pass also values outside of this interval and you
     * can Get result (extrapolation?)
     */
    SIMD_INLINE IVector4D<T> Lerp(T fact, const IVector4D<T>& r) const
    {
        return (*this) + (r - (*this)) * fact;
    }





    //! Returns the angle (in radians) between the two (normalized or unnormalized) vectors 'lhs' and 'rhs'.
    SIMD_INLINE T AngleBetween( const IVector4D<T> &rhs ) const
    {
        IVector4D<T> lhs(*this);
        T dotProduct = lhs.Dot(rhs);
        T vectorsMagnitude = (lhs.Length()) * (rhs.Length());
        T angle = IACos(dotProduct / vectorsMagnitude);
        if( is_nan(angle)) return 0;
        return (angle);
    }

    //-------------[ conversion ]-----------------------------

    /**
     * Conversion to pointer operator
     * @return Pointer to internally stored (in management of class IVector4D<T>)
     * used for passing IVector4D<T> values to gl*4[fd] functions.
     */
    SIMD_INLINE operator T*()
    {
        return &x;
    }

    /**
     * Conversion to pointer operator
     * @return Constant Pointer to internally stored (in management of class IVector4D<T>)
     * used for passing IVector4D<T> values to gl*4[fd] functions.
     */
    SIMD_INLINE operator const T*() const
    {
        return &x;
    }

    //-------------[ output operator ]------------------------
    /**
    * Output to stream operator
    * @param lhs Left hand side argument of operator (commonly ostream instance).
    * @param rhs Right hand side argument of operator.
    * @return Left hand side argument - the ostream object passed to operator.
    */
    friend std::ostream& operator<<(std::ostream& lhs, const IVector4D<T>& rhs)
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


public:
    /**
     * The multiplicitive identity vector
     */
    static const IVector4D<T> IDENTITY;

    /**
     * The additive identity vector.
     */
    static const IVector4D<T> ZERO;
    /**
     * The identity vector X.
     */
    static const IVector4D<T> X;
    /**
     * The identity vector Y.
     */
    static const IVector4D<T> Y;

    /**
     * The identity vector Z.
     */
    static const IVector4D<T> Z;

    /**
     * The identity vector W.
     */
    static const IVector4D<T> W;

};


template<class T> const IVector4D<T> IVector4D<T>::IDENTITY(1.0, 1.0, 1.0, 1.0);
template<class T> const IVector4D<T> IVector4D<T>::ZERO(0.0, 0.0, 0.0, 0.0);
template<class T> const IVector4D<T> IVector4D<T>::X(1.0, 0.0, 0.0, 0.0);
template<class T> const IVector4D<T> IVector4D<T>::Y(0.0, 1.0, 0.0, 0.0);
template<class T> const IVector4D<T> IVector4D<T>::Z(0.0, 0.0, 1.0, 0.0);
template<class T> const IVector4D<T> IVector4D<T>::W(0.0, 0.0, 0.0, 1.0);


template<class T> const
static IVector4D<T> Cross(const IVector4D<T>& a, const IVector4D<T>& b , const IVector4D<T>& c)
{
    return a.Cross(b,c);
}

template<class T> const
static T Dot(const IVector4D<T>& a, const IVector4D<T>& b)
{
    return a.Dot(b);
}




//--------------------------------------
// Typedef shortcuts for 4D vector
//-------------------------------------

using IVector4r    = IVector4D<Real>;
using IVector4f    = IVector4D<float>;
using IVector4d    = IVector4D<double>;
using IVector4i    = IVector4D<std::int32_t>;
using IVector4ui   = IVector4D<std::uint32_t>;
using IVector4b    = IVector4D<std::int8_t>;
using IVector4ub   = IVector4D<std::uint8_t>;


} /* namespace */


#endif // IVECTOR4D_H
