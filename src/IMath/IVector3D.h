/********************************************************************************
 *
 * IVector3D.h
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


#ifndef IVECTOR3D_H_
#define IVECTOR3D_H_


#include "IReal.h"
#include "IFunc.h"
#include "IVector2D.h"


namespace IMath
{


template<class T> class  IMatrix3x3;

/**
 * Class for three dimensional vector.
 * There are four ways of accessing vector components.
 * Let's have <code>Vector3f v</code>, you can either:
 * <ul>
 * 	<li>access as position (x,y,z) &mdash; <code>v.x = v.y = v.z = 1;</code></li>
 * 	<li>access as texture coordinate (s,t,u) &mdash; <code>v.s = v.t = v.u = 1;</code></li>
 * 	<li>access as color (r,g,b) &mdash; <code>v.r = v.g = v.b = 1;</code></li>
 * 	<li>access via operator[] &mdash; <code>v[0] = v[1] = v[2] = 1;</code></li>
 * </ul>
 **/
template<class T>
class  IVector3D
{

public:

    //! Specifies the typename of the scalar components.
    using ScalarType = T;

    //! Specifies the number of vector components.
    static const std::size_t components = 3;

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

        /**
         * First element of vector, alias for R-coordinate.
         * For color notation.
         */
        T r;
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
        /**
         * Second element of vector, alias for G-coordinate.
         * For color notation.
         */
        T g;
    };

    union
    {
        /**
         * Third element of vector, alias for Z-coordinate.
         */
        T z;

        /**
         * Third element of vector, alias for W-coordinate.
         * For textures notation.
         */
        T w;
        /**
         * Third element of vector, alias for B-coordinate.
         * For color notation.
         */
        T b;
    };


public:

    //----------------[ constructors ]--------------------------
    /**
     * Creates and Sets to (0,0,0)
     */
    SIMD_INLINE IVector3D()
        : x(0), y(0), z(0)
    {
    }




    SIMD_INLINE explicit IVector3D(const T& scalar)
        : x ( scalar ) ,
          y ( scalar ) ,
          z ( scalar )
    {

    }


    SIMD_INLINE IVector3D(T* data)
    {
         &x = data;
    }

    /**
     * Creates and Sets to (x,y,z)
     * @param nx initial x-coordinate value
     * @param ny initial y-coordinate value
     * @param nz initial z-coordinate value
     */
    SIMD_INLINE IVector3D(T nx, T ny, T nz)
        : x(nx), y(ny), z(nz)
    {
    }

    /**
     * Copy constructor.
     * @param src Source of data for new created IVector3D instance.
     */
    SIMD_INLINE IVector3D(const IVector3D<T>& src)
        : x(src.x), y(src.y), z(src.z)
    {
    }

    /**
     * Copy casting constructor.
     * @param src Source of data for new created IVector3D instance.
     */
    template<class FromT>
    SIMD_INLINE IVector3D(const IVector3D<FromT>& src)
        : x(static_cast<T>(src.x)),
          y(static_cast<T>(src.y)),
          z(static_cast<T>(src.z))
    {
    }


    //---------------------- Methods ---------------------//

    SIMD_INLINE void SetToZero()
    {
        x = T(0);
        y = T(0);
        z = T(0);
    }

    SIMD_INLINE bool IsZero() const
    {
        return IApproxEqual<T>( LengthSquare(), T(0.0) );
    }


    SIMD_INLINE void SetAllValues(T newX, T newY, T newZ)
    {
        x = newX;
        y = newY;
        z = newZ;
    }


    SIMD_INLINE T GetX() const { return x; }
    SIMD_INLINE T GetY() const { return y; }
    SIMD_INLINE T GetZ() const { return z; }


    SIMD_INLINE void SetX(T _x) { x = _x; }
    SIMD_INLINE void SetY(T _y) { y = _y; }
    SIMD_INLINE void SetZ(T _z) { z = _z; }

    SIMD_INLINE void SetXY( const IVector2D<T>& v)
    {
        x = v.x;
        y = v.y;
    }


    SIMD_INLINE IVector2D<T> GetXY() const { return IVector2D<T>(x,y); }


    //----------------[ access operators ]-------------------
    /**
     * Copy operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T> operator=(const IVector3D<T>& rhs)
    {
        x = rhs.x;
        y = rhs.y;
        z = rhs.z;
        return *this;
    }

    /**
     * Copy casting operator.
     * @param rhs Right hand side argument of binary operator.
     */
    template<class FromT>
    SIMD_INLINE IVector3D<T> operator=(const IVector3D<FromT>& rhs)
    {
        x = static_cast<T>(rhs.x);
        y = static_cast<T>(rhs.y);
        z = static_cast<T>(rhs.z);
        return *this;
    }

    /**
     * Array access operator
     * @param n Array index
     * @return For n = 0, reference to x coordinate, n = 1
     * reference to y, else reference to z
     * y coordinate.
     */
    SIMD_INLINE T & operator[](int n)
    {
        static_assert(sizeof(*this) == sizeof(T[components]), "");
        assert(n >= 0 && n < 3);
        return (&x)[n];
    }

    /**
     * Constant array access operator
     * @param n Array index
     * @return For n = 0, reference to x coordinate, n = 1
     * reference to y, else reference to z
     * y coordinate.
     */
    SIMD_INLINE const T & operator[](int n) const
    {
        //static_assert(sizeof(*this) == sizeof(T[components]), "");
        assert(n >= 0 && n < 3);
        return (&x)[n];
    }

    //---------------[ vector arithmetic operator ]--------------
    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T> operator+(const IVector3D<T>& rhs) const
    {
        return IVector3D<T>(x + rhs.x, y + rhs.y, z + rhs.z);
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T> operator-(const IVector3D<T>& rhs) const
    {
        return IVector3D<T>(x - rhs.x, y - rhs.y, z - rhs.z);
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T> operator*(const IVector3D<T>& rhs) const
    {
        return IVector3D<T>(x * rhs.x, y * rhs.y, z * rhs.z);
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T> operator/(const IVector3D<T>& rhs) const
    {
        return IVector3D<T>(x / rhs.x, y / rhs.y, z / rhs.z);
    }

    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T>& operator+=(const IVector3D<T>& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T>& operator-=(const IVector3D<T>& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T>& operator*=(const IVector3D<T>& rhs)
    {
        x *= rhs.x;
        y *= rhs.y;
        z *= rhs.z;
        return *this;
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T>& operator/=(const IVector3D<T>& rhs)
    {
        x /= rhs.x;
        y /= rhs.y;
        z /= rhs.z;
        return *this;
    }


    //------------------------------ Dot . Cross ----------------------------------------//

    /**
     * Dot product of two vectors.
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE T Dot(const IVector3D<T>& rhs) const
    {
        return x * rhs.x + y * rhs.y + z * rhs.z;
    }

    /**
     * Cross product operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T> Cross(const IVector3D<T>& rhs) const
    {
        return IVector3D<T>(y * rhs.z - rhs.y * z,
                            z * rhs.x - rhs.z * x,
                            x * rhs.y - rhs.x * y);
    }

    //--------------[ scalar vector operator ]--------------------
    /**
    * Addition operator
    * @param rhs Right hand side argument of binary operator.
    */
    SIMD_INLINE IVector3D<T> operator+(T rhs) const
    {
        return IVector3D<T>(x + rhs, y + rhs, z + rhs);
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T> operator-(T rhs) const
    {
        return IVector3D<T>(x - rhs, y - rhs, z - rhs);
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T> operator*(T rhs) const
    {
        return IVector3D<T>(x * rhs, y * rhs, z * rhs);
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T> operator/(T rhs) const
    {
        return IVector3D<T>(x / rhs, y / rhs, z / rhs);
    }


    //------------------------------ Friends ----------------------------------------//

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    friend SIMD_INLINE  IVector3D<T> operator*(T number, const IVector3D<T>& vector)
    {
        return IVector3D<T>(number * vector.x, number * vector.y, number * vector.z);
    }


    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    friend SIMD_INLINE  IVector3D<T> operator/( T number , const IVector3D<T>& vector )
    {
        return IVector3D<T>(vector.x / number, vector.y / number, vector.z / number);
    }


    //------------------------------------------------------------------------------//

    /**
     * Addition operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T>& operator+=(T rhs)
    {
        x += rhs;
        y += rhs;
        z += rhs;
        return *this;
    }

    /**
     * Subtraction operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T>& operator-=(T rhs)
    {
        x -= rhs;
        y -= rhs;
        z -= rhs;
        return *this;
    }

    /**
     * Multiplication operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T>& operator*=(T rhs)
    {
        x *= rhs;
        y *= rhs;
        z *= rhs;
        return *this;
    }

    /**
     * Division operator
     * @param rhs Right hand side argument of binary operator.
     */
    SIMD_INLINE IVector3D<T>& operator/=(T rhs)
    {
        x /= rhs;
        y /= rhs;
        z /= rhs;
        return *this;
    }

    //--------------[ Equality operator ]------------------------
    /**
     * Equality test operator
     * @param rhs Right hand side argument of binary operator.
     * @note Test of equality is based of threshold EPSILON value. To be two
     * values equal, must satisfy this condition | lhs.x - rhs.y | < EPSILON,
     * same for y-coordinate, and z-coordinate.
     */
    SIMD_INLINE bool operator==(const IVector3D<T>& rhs) const
    {
        return IAbs(x - rhs.x) < MACHINE_EPSILON &&
               IAbs(y - rhs.y) < MACHINE_EPSILON &&
               IAbs(z - rhs.z) < MACHINE_EPSILON;
    }

    /**
     * Inequality test operator
     * @param rhs Right hand side argument of binary operator.
     * @return not (lhs == rhs) :-P
     */
    SIMD_INLINE bool operator!=(const IVector3D<T>& rhs) const
    {
        return !(*this == rhs);
    }

    //-------------[ unary operations ]--------------------------
    /**
     * Unary negate operator
     * @return negated vector
     */
    SIMD_INLINE IVector3D<T> operator-() const
    {
        return IVector3D<T>(-x, -y, -z);
    }

    //-------------[ size operations ]---------------------------
    /**
     * Get length of vector.
     * @return lenght of vector
     */
    SIMD_INLINE T Length() const
    {
        return ISqrt(x * x + y * y + z * z);
    }

    /**
     * Return square of length.
     * @return length ^ 2
     * @note This method is faster then length(). For comparison
     * of length of two vector can be used just this value, instead
     * of more expensive length() method.
     */
    SIMD_INLINE T LengthSquare() const
    {
        return x * x + y * y + z * z;
    }


    SIMD_INLINE int GetMinAxis() const
    {
        return (x < y ? (x < z ? 0 : 2) : (y < z ? 1 : 2));
    }

    SIMD_INLINE int GetMaxAxis() const
    {
        return (x < y ? (y < z ? 2 : 1) : (x < z ? 2 : 0));
    }

    SIMD_INLINE T GetMinValue() const
    {
        return IMin(IMin(x, y), z);
    }

    SIMD_INLINE T GetMaxValue() const
    {
        return IMax(IMax(x, y), z);
    }

    //*********************************************//

    /**
     * Normalize vector
     */
    SIMD_INLINE void Normalize()
    {
        T s = Length();
        x /= s;
        y /= s;
        z /= s;
    }


    /**
     * Normalize unit vector
     */
    SIMD_INLINE IVector3D<T> GetUnit() const
    {
        T lengthVector = Length();
        if (lengthVector < MACHINE_EPSILON)
        {
            return *this;
        }
        // Compute and return the unit vector
        T lengthInv = T(1.0) / lengthVector;
        return IVector3D<T>( x * lengthInv,
                             y * lengthInv,
                             z * lengthInv);
    }



    /**
     * Normalize Unit vector (popular name to methods)
     */
    SIMD_INLINE IVector3D<T> Normalized() const
    {
        T lengthVector = Length();
        if (lengthVector < MACHINE_EPSILON)
        {
            return *this;
        }
        // Compute and return the unit vector
        T lengthInv = T(1.0) / lengthVector;
        return IVector3D<T>( x * lengthInv,
                             y * lengthInv,
                             z * lengthInv);
    }

    /**
     * Inverse vector
     */
     SIMD_INLINE IVector3D<T> Inverse() const
     {
         return IVector3D<T>( T(1.0/x) , T(1.0/y) , T(1.0/z) );
     }


    /**
    * Transpose Vector
    */
    SIMD_INLINE IVector3D<T> Transpose() const
    {
          return IVector3D<T>(z,y,x);
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


    /**
    * Orthogonal unit vector
    */
    SIMD_INLINE IVector3D<T> GetOneUnitOrthogonalVector() const
    {
        assert(Length() > MACHINE_EPSILON);

        // Get the minimum element of the vector
        IVector3D<T> vectorAbs(IAbs(x), IAbs(y), IAbs(z));
        int minElement = vectorAbs.GetMinAxis();

        if (minElement == 0)
        {
            return IVector3D<T>(0.0, -z, y) / ISqrt(y*y + z*z);
        }
        else if (minElement == 1)
        {
            return IVector3D<T>(-z, 0.0, x) / ISqrt(x*x + z*z);
        }
        else
        {
            return IVector3D<T>(-y, x, 0.0) / ISqrt(x*x + y*y);
        }
    }


    //------------[ other operations ]---------------------------
    /**
     * Rotate vector around three axis.
     * @param ax Angle (in degrees) to be rotated around X-axis.
     * @param ay Angle (in degrees) to be rotated around Y-axis.
     * @param az Angle (in degrees) to be rotated around Z-axis.
     */
    SIMD_INLINE void Rotate(T ax, T ay, T az)
    {
        T a = ICos(/*IDegreesToRadians*/(ax));
        T b = ISin(/*IDegreesToRadians*/(ax));
        T c = ICos(/*IDegreesToRadians*/(ay));
        T d = ISin(/*IDegreesToRadians*/(ay));
        T e = ICos(/*IDegreesToRadians*/(az));
        T f = ISin(/*IDegreesToRadians*/(az));
        T nx = c * e * x - c * f * y + d * z;
        T ny = (a * f + b * d * e) * x + (a * e - b * d * f) * y - b * c * z;
        T nz = (b * f - a * d * e) * x + (a * d * f + b * e) * y + a * c * z;
        x = nx;
        y = ny;
        z = nz;

    }



    /**
    \brief Rotates the specified vector 'vec' around the vector 'axis'.
    \param[in] vec Specifies the vector which is to be rotated.
    \param[in] axis Specifies the axis vector to rotate around.
    \param[in] angle Specifies the rotation angle (in radians).
    \return The new rotated vector.
    */
    SIMD_INLINE IVector3D<T> RotateVectorAroundAxis(IVector3D<T> axis, T angle)
    {
        const IVector3D<T>& vec(*this);

        axis.Normalize();

        auto s       = std::sin(angle);
        auto c       = std::cos(angle);
        auto cInv    = T(1) - c;

        IVector3D<T> row0, row1, row2;

        row0.x = axis.x*axis.x + c*(T(1) - axis.x*axis.x);
        row0.y = axis.x*axis.y*cInv - s*axis.z;
        row0.z = axis.x*axis.z*cInv + s*axis.y;

        row1.x = axis.x*axis.y*cInv + s*axis.z;
        row1.y = axis.y*axis.y + c*(T(1) - axis.y*axis.y);
        row1.z = axis.y*axis.z*cInv - s*axis.x;

        row2.x = axis.x*axis.z*cInv - s*axis.y;
        row2.y = axis.y*axis.z*cInv + s*axis.x;
        row2.z = axis.z*axis.z + c*(T(1) - axis.z*axis.z);

        return IVector3D<T>( vec.Dot(row0),
                             vec.Dot(row1),
                             vec.Dot(row2) );
    }

    /**
     * Linear interpolation of two vectors
     * @param fact Factor of interpolation. For translation from positon
     * of this vector to vector r, values of factor goes from 0.0 to 1.0.
     * @param r Second Vector for interpolation
     * @note However values of fact parameter are reasonable only in interval
     * [0.0 , 1.0], you can pass also values outside of this interval and you
     * can Get result (extrapolation?)
     */
    SIMD_INLINE IVector3D<T> Lerp(T fact, const IVector3D<T>& r) const
    {
        return (*this) + (r - (*this)) * fact;
    }



    //! Returns the angle (in radians) between the two (Normalized or unNormalized) vectors 'lhs' and 'rhs'.
    SIMD_INLINE T AngleBetween( const IVector3D<T> &rhs ) const
    {
        IVector3D<T> lhs(*this);
        T dotProduct = lhs.Dot(rhs);
        T vectorsMagnitude = (lhs.Length()) * (rhs.Length());
        T angle = IACos(dotProduct / vectorsMagnitude);
        if( is_nan(angle)) return 0;
        return (angle);
    }




    //-------------[ conversion ]-----------------------------

    /**
     * Conversion to pointer operator
     * @return Pointer to internally stored (in management of class IVector3D<T>)
     * used for passing IVector3D<T> values to gl*3[fd] functions.
     */
    SIMD_INLINE operator T*()
    {
        return &x;
    }

    /**
     * Conversion to pointer operator
     * @return Constant Pointer to internally stored (in management of class IVector3D<T>)
     * used for passing IVector3D<T> values to gl*3[fd] functions.
     */
    SIMD_INLINE operator const T*() const
    {
        return &x;
    }


    SIMD_INLINE int GetHashCode() const
    {
        // Overflow is fine, just wrap
        int hash = (int) 2166136261;
        // Suitable nullity checks etc, of course :)
        hash = (hash * 16777619) ^ std::hash<T>{}(x);
        hash = (hash * 16777619) ^ std::hash<T>{}(y);
        hash = (hash * 16777619) ^ std::hash<T>{}(z);
        return hash;

    }

    //-------------[ output operator ]------------------------
    /**
         * Output to stream operator
         * @param lhs Left hand side argument of operator (commonly ostream instance).
         * @param rhs Right hand side argument of operator.
         * @return Left hand side argument - the ostream object passed to operator.
         */
    friend std::ostream& operator<<(std::ostream& lhs, const IVector3D<T> rhs)
    {
        lhs << "[" << rhs[0] << "," << rhs[1] << "," << rhs[2] << "]";
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



    //========================== plugins =========================//




    /// <summary>
    /// Returns a <see cref="Vector3"/> containing the 3D Cartesian coordinates of a point specified in Barycentric coordinates relative to a 3D triangle.
    /// </summary>
    /// <param name="value1">A <see cref="Vector3"/> containing the 3D Cartesian coordinates of vertex 1 of the triangle.</param>
    /// <param name="value2">A <see cref="Vector3"/> containing the 3D Cartesian coordinates of vertex 2 of the triangle.</param>
    /// <param name="value3">A <see cref="Vector3"/> containing the 3D Cartesian coordinates of vertex 3 of the triangle.</param>
    /// <param name="amount1">Barycentric coordinate b2, which expresses the weighting factor toward vertex 2 (specified in <paramref name="value2"/>).</param>
    /// <param name="amount2">Barycentric coordinate b3, which expresses the weighting factor toward vertex 3 (specified in <paramref name="value3"/>).</param>
    /// <param name="result">When the method completes, contains the 3D Cartesian coordinates of the specified point.</param>
   static IVector3D<T> Barycentric(const IVector3D<T> &value1, const IVector3D<T> &value2, const IVector3D<T>& value3, T amount1, T amount2)
    {
        return IVector3D<T>((value1.x + (amount1 * (value2.x - value1.x))) + (amount2 * (value3.x - value1.x)),
                            (value1.y + (amount1 * (value2.y - value1.y))) + (amount2 * (value3.y - value1.y)),
                            (value1.z + (amount1 * (value2.z - value1.z))) + (amount2 * (value3.z - value1.z)));
    }


    /**
     * Compute the determinant of a matrix whose columns are three given vectors.
     * Useful property: det(a, b, c) = det(c, a, b) = det(b, c, a).
     */
    static SIMD_INLINE T Determinant(const IVector3D<T>& a, const IVector3D<T>& b, const IVector3D<T>& c)
    {
        return a.Dot(b.Cross(c));
    }


    static SIMD_INLINE IVector3D<T> Clamp(const IVector3D<T>& vector, T maxLength)
    {
        if (vector.LengthSquare() > maxLength * maxLength)
        {
            return vector.GetUnit() * maxLength;
        }

        return vector;
    }


    static T SIMD_INLINE AngleSigned(IVector3D<T> v1, IVector3D<T> v2, IVector3D<T> normal)
    {
        return IAtan2( normal.Dot(v1.Cross(v2)), v1.Dot(v2));
    }


    /**
     * @brief triNormal
     * @param V0
     * @param V1
     * @param V2
     * @return Triangle_FaceNormal
     */
    static  SIMD_INLINE  IVector3D<T> triNormal( const IVector3D<T>& V0,
                                                 const IVector3D<T>& V1,
                                                 const IVector3D<T>& V2)
    {
        IVector3D<T> Norm;
        IVector3D<T> E = V1;
        IVector3D<T> F = V2;
        E -= V0;
        F -= V1;
        Norm = E.Cross(F);
        Norm.Normalize();
        return Norm;
    }


    /**
     * @brief triArea
     * @param V0
     * @param V1
     * @param V2
     * @return Triangle_Area
     */
    static  SIMD_INLINE  T Area(const IVector3D<T>& V0,
                                const IVector3D<T>& V1,
                                const IVector3D<T>& V2)
    {
            return T(0.5) * (V1 - V0).Cross(V2 - V0).Length();
    }

    /**
     * @brief BiUnitGrammSchmidt
     * @param n -axis
     * @param p -ortogonal_left
     * @param q -ortogonal_up
     * Gramm schmidt process
     * https://www.math.hmc.edu/calculus/tutorials/gramschmidt/gramschmidt.pdf
     */
    static SIMD_INLINE void BiUnitGrammSchmidt(const IVector3D<T>& n_axis , IVector3D<T>& p, IVector3D<T>& q )
    {
        p = n_axis.GetOneUnitOrthogonalVector();
        q = n_axis.Cross(p);
    }


#define btRecipSqrt(x)  (1.0/ISqrt(x))
#define SIMDSQRT12      (0.7071067811865475244008443621048490)
    /// Bullet physics version (Gramm schmidt process)
    static SIMD_INLINE void BiUnitOrthogonalVector(const IVector3D<T>& n, IVector3D<T>& p, IVector3D<T>& q)
    {
        if (IAbs(n[2]) > SIMDSQRT12)
        {
            // choose p in y-z plane
            T a = n[1]*n[1] + n[2]*n[2];
            T k = btRecipSqrt(a);
            p[0] = T(0);
            p[1] = -n[2]*k;
            p[2] = n[1]*k;
            // Set q = n x p
            q[0] = a*k;
            q[1] = -n[0]*p[2];
            q[2] = n[0]*p[1];
        }
        else
        {
            // choose p in x-y plane
            T a = n[0]*n[0] + n[1]*n[1];
            T k = btRecipSqrt(a);
            p[0] = -n[1]*k;
            p[1] = n[0]*k;
            p[2] = T(0);
            // Set q = n x p
            q[0] = -n[2]*p[1];
            q[1] = n[2]*p[0];
            q[2] = a*k;
        }
    }



public:
    /**
     * The multiplicitive identity vector
     */
    static const IVector3D<T> IDENTITY;
    /**
     * The additive identity vector.
     */
    static const IVector3D<T> ZERO;
    /**
     * The identity vector X.
     */
    static const IVector3D<T> X;
    /**
     * The identity vector Y.
     */
    static const IVector3D<T> Y;
    /**
     * The identity vector Z.
     */
    static const IVector3D<T> Z;


};

template<class T> const IVector3D<T> IVector3D<T>::IDENTITY(1.0, 1.0, 1.0);
template<class T> const IVector3D<T> IVector3D<T>::ZERO(0.0, 0.0, 0.0);
template<class T> const IVector3D<T> IVector3D<T>::X(1.0, 0.0, 0.0);
template<class T> const IVector3D<T> IVector3D<T>::Y(0.0, 1.0, 0.0);
template<class T> const IVector3D<T> IVector3D<T>::Z(0.0, 0.0, 1.0);


template<class T> const
static IVector3D<T> Cross(const IVector3D<T>& a, const IVector3D<T>& b)
{
    return a.Cross(b);
}

template<class T> const
static T Dot(const IVector3D<T>& a, const IVector3D<T>& b)
{
    return a.Dot(b);
}

//--------------------------------------
// Typedef shortcuts for 3D vector
//-------------------------------------

using IVector3r    = IVector3D<Real>;
using IVector3f    = IVector3D<float>;
using IVector3d    = IVector3D<double>;
using IVector3i    = IVector3D<std::int32_t>;
using IVector3ui   = IVector3D<std::uint32_t>;
using IVector3b    = IVector3D<std::int8_t>;
using IVector3ub   = IVector3D<std::uint8_t>;

} /* namespace */



#endif /* IVECTOR3D_H_ */
