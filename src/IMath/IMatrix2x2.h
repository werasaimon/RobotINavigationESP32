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


#ifndef IMATRIX2X2_H
#define IMATRIX2X2_H



#include "IVector2D.h"


namespace IMath
{


template<class T> class IMatrix2x2
{
    public:


        //! Specifies the typename of the scalar components.
        using ScalarType = T;

        //! Specifies the number of vector components.
        static const std::size_t components = 2*2;


    private:
	 //-------------------- Attributes --------------------//

	    union
	    {
            T               mData[4];
            IVector2D<T>    mRows[2];
	    };



    public:
        /**
         * Default constructor.
         *
         * Constructs the identity matrix.
         */
        SIMD_INLINE IMatrix2x2()
        : mData{1.0f, 0.0f,
                0.0f, 1.0f}
        {
                // Nothing to do.
        }


        /// Constructor
        SIMD_INLINE IMatrix2x2(T value)
        {
            setAllValues(value, value,
                         value, value);
        }


        /**
         * Constructor.
         *
         * Constructs the matrix from the passed array.
         *
         * @param arr Array of floating point values in row-major order.
         */
        SIMD_INLINE IMatrix2x2(const T arr[components])
        {
            std::memcpy(mData, arr, components * sizeof(T));
        }


        /**
         * Constructor.
         *
         * Constructs the matrix from the specified entries.
         *
         * @param entry00 Entry at row 0 column 0.
         * @param entry01 Entry at row 0 column 1.
         * @param entry10 Entry at row 1 column 0.
         * @param entry11 Entry at row 1 column 1.
         * @param entry20 Entry at row 2 column 0.
         */
        SIMD_INLINE IMatrix2x2(T entry00, T entry01,
                               T entry10, T entry11)
          : mData{entry00, entry01,
                  entry10, entry11}
        {
                // Nothing to do.
        }

        /**
         * Copy constructor.
         *
         * @param other The other matrix to copy.
         */
        SIMD_INLINE IMatrix2x2(const IMatrix2x2<T>& other) = default;



        /**
        * Copy casting constructor.
        * @param src Data source for new created instance of IMatrix2x2
        */
        template<class FromT>
        SIMD_INLINE IMatrix2x2(const IMatrix2x2<FromT>& src)
        {
            for (int i = 0; i < components; i++)
            {
                mData[i] = static_cast<T>(src.mData[i]);
            }
        }



        //---------------------- Methods ---------------------//

        /// Set the matrix to the identity matrix
        SIMD_INLINE void SetToIdentity()
        {
        	// Initialize all values in the matrix to identity
            SetAllValues(1.0, 0.0,
                         0.0, 1.0);
        }


        /// Set the matrix to zero
        SIMD_INLINE void SetToZero()
        {
            mRows[0].SetToZero();
            mRows[1].SetToZero();
        }


        /// Set all the values in the matrix
        SIMD_INLINE void SetAllValues(T a1, T a2,
                                      T b1, T b2)
        {
        	mRows[0][0] = a1; mRows[0][1] = a2;
            mRows[1][0] = b1; mRows[1][1] = b2;
        }



        /**
        * normalized matrix to be value matrix2x2
        */
        void OrthoNormalize()
        {
           mRows[0].Normalize();
           mRows[1].Normalize();
        }

        /**
        * Normalized matrix to be value matrix2x2
        * Return matrix2x2
        */
        SIMD_INLINE IMatrix2x2<T> OrthoNormalized() const
        {
           IMatrix2x2<T> res(*this);
           res.OrthoNormalize();
           return res;
        }


        /**
        * Matrix entry accessor operator.
        *
        * @note Entry indicies are in the range 0 <= `index` <= 8.
        *
        * @param index Index for the entry to return.
        * @return Entry at position `index`.
        */
        SIMD_INLINE T operator[](std::size_t index) const
        {
            assert(index < components);
            return mData[index];
        }

        /**
         * Equality operator.
         *
         * @param A The first matrix.
         * @param B The second matrix.
         * @return True if the two supplied matrices are equal. False otherwise.
         */
        friend SIMD_INLINE bool operator==(const IMatrix2x2<T>& A, const IMatrix2x2<T>& B)
        {
            const T epsilon = MACHINE_EPSILON;
            for (int i = 0; i < components; ++i)
            {
                if (IAbs(A[i] - B[i]) > epsilon) return false;
            }

            return true;
        }

        /**
         * Non-equality operator.
         *
         * @param A The first matrix.
         * @param B The second matrix.
         * @return True if the two supplied matrices are not equal. False
         * otherwise.
         */
        friend SIMD_INLINE  bool operator!=(const IMatrix2x2<T>& A, const IMatrix2x2<T>& B)
        {
            return !(A == B);
        }



        //---------------------[ assignment operations ]---------------------------------

        /**
        * Copy operator
        * @param rhs Right hand side argument of binary operator.
        */
        SIMD_INLINE IMatrix2x2<T>& operator=(const IMatrix2x2<T>& rhs)
        {
            std::memcpy(mData, rhs.mData, sizeof(T) * components);
            return *this;
        }

        /**
        * Copy casting operator
        * @param rhs Right hand side argument of binary operator.
        */
        template<class FromT>
        SIMD_INLINE IMatrix2x2<T>& operator=(const IMatrix2x2<FromT>& rhs)
        {
            for (int i = 0; i < components; i++)
            {
                mData[i] = static_cast<T>(rhs.mData[i]);
            }
            return *this;
        }

        /**
        * Copy operator
        * @param rhs Right hand side argument of binary operator.
        */
        SIMD_INLINE IMatrix2x2<T>& operator=(const T* rhs)
        {
            std::memcpy(mData, rhs, sizeof(T) * components);
            return *this;
        }


        //-------------[ conversion data ]-----------------------------

        /**
         * Conversion to pointer operator
         * @return Pointer to internally stored (in management of class IMatrix2x2<T>)
         * used for passing IMatrix2x2<T> values to gl*[fd]v functions.
         */
        SIMD_INLINE operator T*()
        {
            return  &mRows[0][0];
        }

        /**
         * Conversion to pointer operator
         * @return Constant Pointer to internally stored (in management of class IMatrix2x2<T>)
         * used for passing IMatrix2x2<T> values to gl*[fd]v functions.
         */
        SIMD_INLINE operator const T*() const
        {
            return  &mRows[0][0];
        }


        /**
         * Conversion to pointer operator
         */
        SIMD_INLINE const T* GetData() const
        {
        	return &mRows[0][0];
        }

        /**
         * Conversion to pointer operator
         */
        SIMD_INLINE T* GetData()
        {
        	return &mRows[0][0];
        }




        /// Overloaded operator to read element of the matrix.
        SIMD_INLINE const IVector2D<T>& operator[](int row) const
        {
        	return mRows[row];
        }

        /// Overloaded operator to read/write element of the matrix.
        SIMD_INLINE IVector2D<T>& operator[](int row)
        {
        	return mRows[row];
        }




        /// Return a column
        SIMD_INLINE IVector2D<T> GetColumn(int i) const
        {
        	assert(i>= 0 && i<2);
            return IVector2D<T> (mRows[0][i], mRows[1][i] );
        }

        /// Return a row
        SIMD_INLINE IVector2D<T> GetRow(int i) const
        {
        	assert(i>= 0 && i<2);
        	return mRows[i];
        }


        /**
         * Matrix addition operator.
         *
         * @note Matrix addition is commutative.
         *
         * @param A The first matrix.
         * @param B The second matrix.
         * @return The matrix equal to the sum of `A` and `B`.
         */
        SIMD_INLINE IMatrix2x2<T> operator+(const IMatrix2x2<T>& rhs) const
        {
            IMatrix2x2<T> ret;
            for (unsigned int i = 0; i < 2; i++)
            {
                ret.mRows[i] = mRows[i] + rhs.mRows[i];
            }
        	return ret;

        }



        /**
         * Matrix subtraction operator.
         *
         * @note Matrix subtraction is not commutative.
         *
         * @param lhs The left hand side matrix.
         * @param rhs The right hand side matrix.
         * @return The matrix equal to the `rhs` matrix subtracted from the `lhs`
         * matrix.
         */
        SIMD_INLINE IMatrix2x2<T> operator-(const IMatrix2x2<T>& rhs) const
        {
            IMatrix2x2<T> ret;
            for (int i = 0; i < 2; i++)
            {
                ret.mRows[i] = mRows[i] - rhs.mRows[i];
            }
        	return ret;
        }



        /**
         * Matrix negation operator.
         *
         * @param A The matrix to negate.
         * @return The additive inverse of the matrix `A`.
         */
        friend SIMD_INLINE IMatrix2x2<T> operator-(const IMatrix2x2<T> &A)
        {
            return IMatrix2x2<T>( -A[0], -A[1],
                                  -A[2], -A[3]
        	);
        }

        /**
         * Scalar multiplication operator.
         *
         * Multiplies each entry of a matrix by a given scalar value.
         *
         * @param A The matrix to be multiplied by the given scalar.
         * @param s The scalar value.
         * @return The matrix `A` multiplied by the scalar `s`.
         */
        friend SIMD_INLINE IMatrix2x2<T> operator*(const IMatrix2x2<T>& matrix, const T nb)
        {
            return IMatrix2x2<T>(matrix.mRows[0][0] * nb, matrix.mRows[0][1] * nb,
                                 matrix.mRows[1][0] * nb, matrix.mRows[1][1] * nb);
        }


        /**
         * Scalar multiplication operator.
         *
         * Multiplies each entry of a matrix by a given scalar value.
         *
         * @param s The scalar value.
         * @param A The matrix to be multiplied by the given scalar.
         * @return The matrix `A` multiplied by the scalar `s`.
         */
        friend SIMD_INLINE IMatrix2x2<T> operator*(const T s, const IMatrix2x2<T>& A)
        {
        	return A * s;
        }

        /// Overloaded operator for inveret multiplication with a number
        friend SIMD_INLINE IMatrix2x2<T>  operator/(T nb, const IMatrix2x2<T>& matrix)
        {
            return IMatrix2x2<T>(matrix.mRows[0][0] / nb, matrix.mRows[0][1] / nb,
                                 matrix.mRows[1][0] / nb, matrix.mRows[1][1] / nb);
        }

        /// Overloaded operator for inveret multiplication with a matrix
        friend SIMD_INLINE IMatrix2x2<T>  operator/(const IMatrix2x2<T>& matrix, T nb)
        {
        	return nb / matrix;
        }

        /**
         * Vector multiplication operator.
         *
         * Multiplies the matrix `rhs` on the left by the row vector `lhs`,
         * returning the resulting vector.
         *
         * @param lhs The matrix.
         * @param rhs The row vector.
         * @return The vector `lhs` multiplied by the matrix `rhs` on the right.
         */
        friend SIMD_INLINE IVector2D<T> operator*(const IVector2D<T>& rhs, const IMatrix2x2<T>& lhs)
        {
            T fX = rhs.x;
            T fY = rhs.y;

            IVector2D<T> Point;
            Point.x = ( fX * lhs.mRows[0][0] + fY * lhs.mRows[1][0] );
            Point.y = ( fX * lhs.mRows[0][1] + fY * lhs.mRows[1][1] );

            return Point;
        }

        /**
         * Vector multiplication operator.
         *
         * Multiplies the column vector `rhs` on the left by the matrix `lhs`,
         * returning the resulting vector.
         *
         * @param lhs The matrix.
         * @param rhs The column vector.
         * @return The vector `rhs` multiplied by the matrix `lhs` on the left.
         */
        SIMD_INLINE IVector2D<T> operator*(const IVector2D<T>& rhs) const
        {
        	T fX = rhs.x;
        	T fY = rhs.y;

            IVector2D<T> Point;
        	Point.x = ( fX * mRows[0][0] + fY * mRows[0][1] );
        	Point.y = ( fX * mRows[1][0] + fY * mRows[1][1] );

        	return Point;

        }



        /**
         * Matrix multiplication operator.
         *
         * @note Matrix multiplication is not commutative.
         *
         * @param lhs The left hand side matrix.
         * @param rhs The right hand side matrix.
         * @return The matrix equal to the product `lhs` x `rhs`.
         */
        SIMD_INLINE IMatrix2x2<T> operator*(IMatrix2x2<T> rhs) const
        {
            IMatrix2x2<T> w;
        	for(int i = 0; i < 2; i++)
        	{
        		for (int j = 0; j < 2; j++)
        		{
        			T n = 0;
        			for (int k = 0; k < 2; k++)
        			{
        				n += rhs.mRows[i][k] * mRows[k][j];
        			}
        			w.mRows[i][j] = n;
        		}
        	}
        	return w;
        }


        //---------------------------[ misc operations ]----------------------------


        /**
         * Returns the determinant of the matrix.
         *
         * @note A square matrix is invertable if and only if its determinant is
         * nonzero.
         *
         * @return The determinant.
         */
        SIMD_INLINE T Determinant() const
        {
        	return mData[0] * mData[3] - mData[1] * mData[2];
        }

        /**
         * Returns a copy of the multiplicitive inverse of this matrix.
         *
         * @return The multiplicitive inverse of this matrix.
         */
        SIMD_INLINE IMatrix2x2<T> Inverse() const
        {
        	// Ensure that the matrix is not singular.
            const T det = Determinant();
        	assert(det != 0.0f);

        	// Return a copy of the inverse of this matrix.
        	const T invDet = 1.0f / det;
            return IMatrix2x2<T>( mData[3] * invDet, -mData[1] * invDet,
                                 -mData[2] * invDet,  mData[0] * invDet );
        }

        void LoadInverse()
        {
            *this = *this->Inverse();
        }

        /**
         * Returns a copy of this matrix transposed so that the rows now form
         * columns.
         *
         * @return Transposed copy of this matrix.
         */
        SIMD_INLINE IMatrix2x2<T> Transpose() const
        {
            IMatrix2x2<T> ret;
        	for (int i = 0; i < 2; i++)
        	{
        		for (int j = 0; j < 2; j++)
        		{
        			ret.mRows[i][j] = mRows[j][i];
        		}
        	}
        	return ret;
        }


        /**
         * Return the matrix with absolute values
         */
        SIMD_INLINE IMatrix2x2<T> AbsoluteMatrix() const
        {
            return IMatrix2x2<T>(IAbs(mRows[0][0]), IAbs(mRows[0][1]),
                                 IAbs(mRows[1][0]), IAbs(mRows[1][1]));
        }


        /**
         * Return the trace of the matrix
         */
        SIMD_INLINE T Trace() const
        {
        	// Compute and return the trace
        	return (mRows[0][0] + mRows[1][1]);
        }


        //--------------------------------------------------------------------//


        /// Overloaded operator for addition with assignment
        SIMD_INLINE IMatrix2x2<T>& operator+=(const IMatrix2x2<T>& matrix)
        {
        	mRows[0] += matrix.mRows[0];
        	mRows[1] += matrix.mRows[1];
        	return *this;
        }

        /// Overloaded operator for substraction with assignment
        SIMD_INLINE IMatrix2x2<T>& operator-=(const IMatrix2x2<T>& matrix)
        {
        	mRows[0] -= matrix.mRows[0];
        	mRows[1] -= matrix.mRows[1];
        	return *this;
        }

        /// Overloaded operator for multiplication with a number with assignment
        SIMD_INLINE IMatrix2x2<T>& operator*=(T nb)
        {
        	mRows[0] *= nb;
        	mRows[1] *= nb;
        	return *this;
        }

        /// Overloaded operator for invert multiplication with a number with assignment
        SIMD_INLINE IMatrix2x2<T> &operator/=(T nb)
        {
        	mRows[0] /= nb;
        	mRows[1] /= nb;
        	return *this;
        }



        //-------------------------- Plugins ----------------------------//


        /// Return a skew-symmetric matrix using a given vector that can be used
        /// to compute Cross product with another vector using matrix multiplication
        static SIMD_INLINE IMatrix2x2<T> ComputeSkewSymmetricMatrixForCrossProduct(const IVector2D<T>& vector)
        {
            return IMatrix2x2<T>(0, -vector.x, vector.y ,0 );
        }




         /**
          * Returns a scaling matrix that scales by `factor` uniformly.
          *
          * @param scale Uniform scale factor.
          * @return Scaling matrix.
          */
         static SIMD_INLINE IMatrix2x2<T> CreateScale(const T _scale)
         {
             static IMatrix2x2 res;

             T _x = T(_scale);
             T _y = T(_scale);

             res.mRows[0] = IVector2D<T>( _x ,  0.f );
             res.mRows[1] = IVector2D<T>( 0.f , _y  );

             return res;
         }



         /**
         * Returns a scaling matrix that scales by `scaleFactors.x` and
         * 'scaleFactors.y' in the x and y axes respectively.
         *
         * @param scaleFactors Scale factors.
         * @return Scaling matrix.
         */
         static SIMD_INLINE IMatrix2x2<T> CreateScale(  T _x , T _y  )
         {
             static IMatrix2x2 res;

             res.mRows[0] = IVector2D<T>( _x  ,  0.f );
             res.mRows[1] = IVector2D<T>( 0.f , _y   );

             return res;
         }



        /**
         * Returns a scaling matrix that scales by `scaleFactors.x` and
         * 'scaleFactors.y' in the x and y axes respectively.
         *
         * @param scaleFactors Scale factors.
         * @return Scaling matrix.
         */
        static SIMD_INLINE IMatrix2x2<T> CreateScale(const IVector2D<T>& _scale )
        {
            static IMatrix2x2 res;

            T _x = T(_scale.x);
            T _y = T(_scale.y);

            res.mRows[0] = IVector2D<T>( _x  ,  0.f );
            res.mRows[1] = IVector2D<T>( 0.f , _y   );

            return res;
        }



        /**
         * Returns a scaling around axis matrix that scales
         * @return axis to scaling matrix.
         */
        static SIMD_INLINE IMatrix2x2<T>  CreateScaleAroundAxis( const IVector2D<T> _axis , T _scale )
        {
            static IMatrix2x2<T> M;
            T bgamma = (_scale - 1.0);

            M[0][0] = 1.0+((bgamma)*((_axis.x * _axis.x)));
            M[1][0] =     ((bgamma)*((_axis.y * _axis.x)));

            M[0][1] =     ((bgamma)*((_axis.x * _axis.y)));
            M[1][1] = 1.0+((bgamma)*((_axis.y * _axis.y)));

            return M;
        }

        /**
         * Returns a rotation matrix that rotates by `angle` radians.
         *
         * @param angle Angle (in radians) for the rotation.
         * @return Rotation matrix that rotates `angle` radians
         * counter-clockwise.
         */
        static SIMD_INLINE IMatrix2x2<T> CreateRotation(const T angle)
        {

            const T cosTheta = ICos(angle);
            const T sinTheta = ISin(angle);
            return IMatrix2x2<T>( cosTheta, -sinTheta,
                                  sinTheta,  cosTheta);
        }

        /**
         * Returns a rotation matrix that represents the sortest rotation from
         * the `fromDirection` to the `toDirection`.
         *
         * @param fromDirection Direction for the matrix to rotate from.
         * @param toDirection Direction for the matrix to rotate to.
         * @return Rotation matrix corresponding to the rotation from
         * `fromDirection` to `toDirection`.
         */
        static SIMD_INLINE IMatrix2x2<T> CreateFromToRotation(const IVector2D<T>& fromDirection, const IVector2D<T>& toDirection)
        {
            assert(fromDirection.LengthSquare() > 0.0f && toDirection.LengthSquare() > 0.0f);

            // Compute the angle between the two vectors.
            const T theta = IVector2D<T>::GetAngleBetween(toDirection);

            // Return the rotation matrix.
            return CreateRotation(theta);
        }



        //----------[ output operator ]----------------------------
        /**
        * Output to stream operator
        * @param lhs Left hand side argument of operator (commonly ostream instance).
        * @param rhs Right hand side argument of operator.
        * @return Left hand side argument - the ostream object passed to operator.
        */
        friend std::ostream& operator <<(std::ostream& lhs, const IMatrix2x2<T>& rhs)
        {
            for (int i = 0; i < 2; i++)
            {
                lhs << "|\t";
                for (int j = 0; j < 2; j++)
                {
                    lhs << rhs[i][j]  << "\t";
                }
                lhs << "|" << std::endl;
            }
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
         * The multiplicitive identity matrix.
         */
        static const IMatrix2x2<T> IDENTITY;

        /**
         * The additive identity matrix.
         */
        static const IMatrix2x2<T> ZERO;

};

template<class T> const IMatrix2x2<T> IMatrix2x2<T>::IDENTITY = IMatrix2x2<T>(1.0, 0.0,0.0, 1.0);

template<class T> const IMatrix2x2<T> IMatrix2x2<T>::ZERO = IMatrix2x2<T>(0.0, 0.0, 0.0, 0.0);


template<class T>  SIMD_INLINE IMatrix2x2<T> operator ^ (const IVector2D<T> lhs , const IVector2D<T> rhs )
{
    return IMatrix2x2<T>( lhs.x * rhs.x , lhs.x * rhs.y,
                          lhs.y * rhs.x , lhs.y * rhs.y);
}


//--------------------------------------
// Typedef shortcuts for Matrix2x2
//-------------------------------------

using IMatrix2x2r    = IMatrix2x2<Real>;
using IMatrix2x2f    = IMatrix2x2<float>;
using IMatrix2x2d    = IMatrix2x2<double>;
using IMatrix2x2i    = IMatrix2x2<std::int32_t>;
using IMatrix2x2ui   = IMatrix2x2<std::uint32_t>;
using IMatrix2x2b    = IMatrix2x2<std::int8_t>;
using IMatrix2x2ub   = IMatrix2x2<std::uint8_t>;


} /* namespace */

#endif // IMATRIX2X2_H
