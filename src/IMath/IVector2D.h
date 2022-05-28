/********************************************************************************
 *
 * IMatrix2x2.h
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

#ifndef IVECTOR2_H
#define IVECTOR2_H

#include "IReal.h"
#include "IFunc.h"


namespace IMath
{

template<class T> class  IMatrix2x2;

/**
 * Class for two dimensional vector.
 * There are three ways of accessing vector components.
 * Let's have <code>Vector2f v</code>, you can either:
 * <ul>
 * 	<li>access as position(x,y) &mdash; <code>v.x = v.y = 3;</code></li>
 * 	<li>access as texture coordinate (s,t) &mdash; <code>v.s = v.t = 3;</code></li>
 * 	<li>access via operator[] &mdash; <code>v[0] = v[1] = 3;</code></li>
 * </ul>
 */


template<class T>
class IVector2D
{
public:


    //! Specifies the typename of the scalar components.
    using ScalarType = T;

    //! Specifies the number of vector components.
    static const std::size_t components = 2;

    //-------------------- Attributes --------------------//

    union
    {
        /**
         * First element of vector, alias for X-coordinate.
         */
        T x;

        /**
         * First element of vector, alias for U-coordinate.
         * For textures notation.
         */
        T u;
    };

    union
    {
        /**
         * Second element of vector, alias for Y-coordinate.
         */
        T y;

        /**
         * Second element of vector, alias for V-coordinate.
         * For textures notation.
         */
        T v;
    };

public:

    //----------------[ constructors ]--------------------------
    /**
     * Creates and Sets to (0,0)
     */
    SIMD_INLINE IVector2D()
        : x(0), y(0)
    {
    }





    SIMD_INLINE explicit IVector2D(const T& scalar)
        : x ( scalar ), y ( scalar )
    {
    }

    /**
     * Creates and Sets to (x,y)
     * @param nx initial x-coordinate value
     * @param ny initial y-coordinate value
     */
    SIMD_INLINE IVector2D(T nx, T ny)
        : x(nx), y(ny)
    {
    }

    //    /**
    //     * Copy constructor.
    //     * @param src Source of data for new created instance.
    //     */
    //    IVector2D(const IVector2D<T>& src)
    //            : x(src.x), y(src.y)
    //    {
    //    }

    /**
     * Copy casting constructor.
     * @param src Source of data for new created instance.
     */
    template<class FromT>
    SIMD_INLINE IVector2D(const IVector2D<FromT>& src)
        : x(static_cast<T>(src.x)),
          y(static_cast<T>(src.y))
    {
    }





    //---------------------- Methods ---------------------//

    SIMD_INLINE void SetToZero()
    {
        x = T(0);
        y = T(0);
    }


    SIMD_INLINE void SetAllValues(T newX, T newY )
    {
        x = newX;
        y = newY;
    }

    SIMD_INLINE T GetX() const { return x; }
    SIMD_INLINE T GetY() const { return y; }

    SIMD_INLINE void SetX(T _x) { x = _x; }
    SIMD_INLINE void SetY(T _y) { y = _y; }


    //----------------[ access operators ]-------------------
    /**
     * Copy casting operator
     * @param rhs Right hand side argument of binary operator.
     */
    template<class FromT>
    SIMD_INLINE IVector2D<T>& operator=(const IVector2D<FromT>& rhs)
    {
        x = static_cast<T>(rhs.x);
        y = static_cast<T>(rhs.y);
        return *this;
    }

    /**
     * Copy operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector2D<T>& operator=(const IVector2D<T>& rhs)
    {
        x = rhs.x;
        y = rhs.y;
        return *this;
    }

    /**
     * Array access operator
     * @param n Array index
     * @return For n = 0, reference to x coordinate, else reference to y
     * y coordinate.
     */
    SIMD_INLINE T& operator[](int n)
    {
        static_assert(sizeof(*this) == sizeof(T[components]), "");
        assert(n >= 0 && n < 2);
        return (&x)[n];
    }

    /**
     * Constant array access operator
     * @param n Array index
     * @return For n = 0, reference to x coordinate, else reference to y
     * y coordinate.
     */
    SIMD_INLINE const T& operator[](int n) const
    {
        static_assert(sizeof(*this) == sizeof(T[components]), "");
        assert(n >= 0 && n < 2);
        return (&x)[n];
    }

    //---------------[ vector arithmetic operator ]--------------
    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector2D<T> operator+(const IVector2D<T>& rhs) const
    {
        return IVector2D<T>(x + rhs.x, y + rhs.y);
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector2D<T> operator-(const IVector2D<T>& rhs) const
    {
        return IVector2D<T>(x - rhs.x, y - rhs.y);
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector2D<T> operator*(const IVector2D<T>& rhs) const
    {
        return IVector2D<T>(x * rhs.x, y * rhs.y);
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector2D<T> operator/(const IVector2D<T>& rhs) const
    {
        return IVector2D<T>(x / rhs.x, y / rhs.y);
    }

    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector2D<T>& operator+=(const IVector2D<T>& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }

    /**
     * Substraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector2D<T>& operator-=(const IVector2D<T>& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector2D<T>& operator*=(const IVector2D<T>& rhs)
    {
        x *= rhs.x;
        y *= rhs.y;
        return *this;
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector2D<T>& operator/=(const IVector2D<T>& rhs)
    {
        x /= rhs.x;
        y /= rhs.y;
        return *this;
    }

    //--------------[ scalar vector operator ]--------------------
    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector2D<T> operator+(T rhs) const
    {
        return IVector2D<T>(x + rhs, y + rhs);
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector2D<T> operator-(T rhs) const
    {
        return IVector2D<T>(x - rhs, y - rhs);
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector2D<T> operator*(T rhs) const
    {
        return IVector2D<T>(x * rhs, y * rhs);
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector2D<T> operator/(T rhs) const
    {
        return IVector2D<T>(x / rhs, y / rhs);
    }

    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector2D<T>& operator+=(T rhs)
    {
        x += rhs;
        y += rhs;
        return *this;
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector2D<T>& operator-=(T rhs)
    {
        x -= rhs;
        y -= rhs;
        return *this;
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector2D<T>& operator*=(T rhs)
    {
        x *= rhs;
        y *= rhs;
        return *this;
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector2D<T>& operator/=(T rhs)
    {
        x /= rhs;
        y /= rhs;
        return *this;
    }


    //------------------------------ Friends ----------------------------------------//

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    friend SIMD_INLINE  IVector2D<T> operator*(T number, const IVector2D<T>& vector)
    {
        return IVector2D<T>(number * vector.x, number * vector.y);
    }


    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    friend SIMD_INLINE  IVector2D<T> operator/( T number , const IVector2D<T>& vector )
    {
        return IVector2D<T>(vector.x / number, vector.y / number);
    }


    //------------------------------ Dot . Cross ----------------------------------------//

    /**
     * Dot product of two vectors.
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE T Dot(const IVector2D<T>& rhs) const
    {
        return x * rhs.x + y * rhs.y;
    }

    /**
     * Cross product operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE T Cross(const IVector2D<T>& rhs) const
    {
        // just calculate the z-component
        return x * rhs.y - y * rhs.x;
    }

    //--------------[ equality operator ]------------------------
    /**
     * Equality test operator
     * @param rhs Right hand side argument of binary operator.
     * @note Test of equality is based of threshold MACHINE_EPSILON value. To be two
     * values equal, must satisfy this condition | lhs.x - rhs.y | < MACHINE_EPSILON,
     * same for y-coordinate.
     */
    SIMD_INLINE bool operator==(const IVector2D<T>& rhs) const
    {
        return (IAbs(x - rhs.x) < MACHINE_EPSILON) &&
               (IAbs(y - rhs.y) < MACHINE_EPSILON);
    }

    /**
     * Inequality test operator
     * @param rhs Right hand side argument of binary operator.
     * @return not (lhs == rhs) :-P
     */
    SIMD_INLINE bool operator!=(const IVector2D<T>& rhs) const
    {
        return !(*this == rhs);
    }

    //-------------[ unary operations ]--------------------------
    /**
     * Unary negate operator
     * @return negated vector
     */
    SIMD_INLINE IVector2D<T> operator-() const
    {
        return IVector2D<T>(-x, -y);
    }

    //-------------[ size operations ]---------------------------

    /**
     * Return square of length.
     * @return length ^ 2
     * @note This method is faster then length(). For comparison
     * of length of two vector can be used just this value, instead
     * of more expensive length()^ 2 method.
     */
    SIMD_INLINE T LengthSquare() const
    {
        return x * x + y * y;
    }


    /**
     * Get length of vector.
     * @return lenght of vector
     */
    SIMD_INLINE T Length() const
    {
        return ISqrt(x * x + y * y);
    }

    /**
     * Normalize vector
     */
    SIMD_INLINE void Normalize()
    {
        T s = Length();
        x /= s;
        y /= s;
    }


    /**
     * Normalize unit vector
     */
    SIMD_INLINE IVector2D<T> GetUnit() const
    {
        T lengthVector = Length();
        if (lengthVector < MACHINE_EPSILON)
        {
            return *this;
        }
        // Compute and return the unit vector
        T lengthInv = T(1.0) / lengthVector;
        return IVector2D<T>( x * lengthInv,
                             y * lengthInv);
    }

    /**
     * Normalize Unit vector (popular name to methods)
     */
    SIMD_INLINE IVector2D<T> Normalized() const
    {
        T lengthVector = Length();
        if (lengthVector < MACHINE_EPSILON)
        {
            return *this;
        }
        // Compute and return the unit vector
        T lengthInv = T(1.0) / lengthVector;
        return IVector2D<T>( x * lengthInv,
                             y * lengthInv);
    }


    /**
     * @Inverse vector
     */
    SIMD_INLINE IVector2D<T> Inverse() const
    {
        return IVector2D<T>( T(1.0/x) , T(1.0/y) );
    }



    /**
    * @Transpose Vector
    */
    SIMD_INLINE IVector2D<T> Transpose() const
    {
          return IVector2D<T>(y,x);
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
    SIMD_INLINE IVector2D<T> Lerp(T fact, const IVector2D<T>& r) const
    {
        return (*this) + (r - (*this)) * fact;
    }



    //! Returns the angle (in radians) between the two (normalized or unnormalized) vectors 'lhs' and 'rhs'.
    SIMD_INLINE T AngleBetween( const IVector2D<T> &rhs ) const
    {
        IVector2D<T> lhs(*this);
        T dotProduct = lhs.Dot(rhs);
        T vectorsMagnitude = (lhs.Length()) * (rhs.Length());
        T angle = IACos(dotProduct / vectorsMagnitude);
        if( is_nan(angle)) return 0;
        return (angle);
    }

    //-------------[ conversion ]-----------------------------
    /**
     * Conversion to pointer operator
     * @return Pointer to internally stored (in management of class IVector2D<T>)
     * used for passing IVector2D<T> values to gl*2[fd] functions.
     */
    SIMD_INLINE operator T*()
    {
        return &x;
    }
    /**
     * Conversion to pointer operator
     * @return Constant Pointer to internally stored (in management of class IVector2D<T>)
     * used for passing IVector2D<T> values to gl*2[fd] functions.
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
    friend std::ostream& operator<<(std::ostream& lhs, const IVector2D<T>& rhs)
    {
        lhs << "[" << rhs[0] << "," << rhs[1] << "]";
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
    static const IVector2D<T> IDENTITY;
    /**
     * The additive identity vector.
     */
    static const IVector2D<T> ZERO;
    /**
     * The identity vector X.
     */
    static const IVector2D<T> X;

    /**
     * The identity vector Y.
     */
    static const IVector2D<T> Y;

};

template<class T> const IVector2D<T> IVector2D<T>::IDENTITY(1.0, 1.0);
template<class T> const IVector2D<T> IVector2D<T>::ZERO(0.0, 0.0);
template<class T> const IVector2D<T> IVector2D<T>::X(1.0, 0.0);
template<class T> const IVector2D<T> IVector2D<T>::Y(0.0, 1.0);


template<class T> const
static T Cross(const IVector2D<T>& a, const IVector2D<T>& b)
{
    return a.Cross(b);
}

template<class T> const
static T Dot(const IVector2D<T>& a, const IVector2D<T>& b)
{
    return a.Dot(b);
}

//--------------------------------------
// Typedef shortcuts for 2D vector
//-------------------------------------

using IVector2r    = IVector2D<Real>;
using IVector2f    = IVector2D<float>;
using IVector2d    = IVector2D<double>;
using IVector2i    = IVector2D<std::int32_t>;
using IVector2ui   = IVector2D<std::uint32_t>;
using IVector2b    = IVector2D<std::int8_t>;
using IVector2ub   = IVector2D<std::uint8_t>;

} /* namespace */



#endif // IVECTOR2_H
