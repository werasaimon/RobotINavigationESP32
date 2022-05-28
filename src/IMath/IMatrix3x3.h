 /********************************************************************************
 *
 * IMatrix3x3.h
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

#ifndef IMATRIX3X3_H_
#define IMATRIX3X3_H_


// Libraries
#include <cassert>
#include "IVector2D.h"
#include "IVector4D.h"



namespace IMath
{



template<class T> class IQuaternion;

//**
// * Class for matrix 3x3.
// * @note Data stored in this matrix are in column major order. This arrangement suits OpenGL.
// * If you're using row major matrix, consider using fromRowMajorArray as way for construction
// * Matrix3<T> instance.
// */
template<class T> class IMatrix3x3
{

    public:


        //! Specifies the typename of the scalar components.
        using ScalarType = T;

        //! Specifies the number of vector components.
        static const std::size_t components = 3*3;

   private:

   //-------------------- Attributes --------------------//

    union
    {
        T             mData[9];
        IVector3D<T>  mRows[3];
    };


    public:


     //--------------------------[ constructors ]-------------------------------

      // Constructor of the class Matrix3x3
     SIMD_INLINE IMatrix3x3()
     {
          // Initialize all values in the matrix to identity
          SetAllValues(1.0, 0.0, 0.0,
                       0.0, 1.0, 0.0,
                       0.0, 0.0, 1.0);
     }

      // Constructor
     SIMD_INLINE IMatrix3x3(T value)
      {
          SetAllValues(value, value, value,
                       value, value, value,
                       value, value, value);
      }

      /**
      * Copy casting constructor.
      * @param src Data source for new created instance of IMatrix3x3
      */
      SIMD_INLINE IMatrix3x3(T a1, T a2, T a3,
                             T b1, T b2, T b3,
                             T c1, T c2, T c3)
      {
                  SetAllValues(a1, a2, a3,
                               b1, b2, b3,
                               c1, c2, c3);
      }


      // Constructor with arguments
      SIMD_INLINE IMatrix3x3( T data[3][3] )
      {
          SetAllValues(data[0][0], data[0][1], data[0][2],
                       data[1][0], data[1][1], data[1][2],
                       data[2][0], data[2][1], data[2][2]);
      }


      // Copy-constructor
      SIMD_INLINE IMatrix3x3(const IMatrix3x3<T>& matrix)
      {
          SetAllValues(matrix.mRows[0][0], matrix.mRows[0][1], matrix.mRows[0][2],
                       matrix.mRows[1][0], matrix.mRows[1][1], matrix.mRows[1][2],
                       matrix.mRows[2][0], matrix.mRows[2][1], matrix.mRows[2][2]);
      }


      /**
      * Copy matrix values from array (these data must be in column
      * major order!)
      */
      SIMD_INLINE IMatrix3x3(const T * dt)
      {
          std::memcpy(mData, dt, sizeof(T) * components);
      }


      /**
      * Copy casting constructor.
      * @param src Data source for new created instance of IMatrix3x3
      */
      template<class FromT>
      SIMD_INLINE IMatrix3x3(const IMatrix3x3<FromT>& src)
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
          SetAllValues(1.0, 0.0, 0.0,
                       0.0, 1.0, 0.0,
                       0.0, 0.0, 1.0);
      }


      /// Set the matrix to zero
      SIMD_INLINE void SetToZero()
      {
          mRows[0].SetToZero();
          mRows[1].SetToZero();
          mRows[2].SetToZero();
      }


      /// Set all the values in the matrix
      SIMD_INLINE void SetAllValues(T a1, T a2, T a3,
                                    T b1, T b2, T b3,
                                    T c1, T c2, T c3)
      {
          mRows[0][0] = a1; mRows[0][1] = a2; mRows[0][2] = a3;
          mRows[1][0] = b1; mRows[1][1] = b2; mRows[1][2] = b3;
          mRows[2][0] = c1; mRows[2][1] = c2; mRows[2][2] = c3;
      }



      /**
       * normalized matrix to be value matrix3x3
       */
      SIMD_INLINE void OrthoNormalize()
      {
         mRows[0].Normalize();
         mRows[1].Normalize();
         mRows[2].Normalize();

//         IVector3D<T> col[3];
//         col[0] = GetColumn(0);
//         col[1] = GetColumn(1);
//         col[2] = GetColumn(2);
//         col[0].Normalize();
//         col[1].Normalize();
//         col[2].Normalize();

//         *this = IMatrix3x3<T>( col[0].x , col[1].x , col[2].x,
//                                col[0].y , col[1].y , col[2].y,
//                                col[0].z , col[1].z , col[2].z);
      }


      /**
      * Normalized matrix to be value matrix3x3
      * Return matrix3x3
      */
      SIMD_INLINE IMatrix3x3<T> OrthoNormalized() const
      {
         IMatrix3x3<T> res(*this);
         res.OrthoNormalize();
         return res;
      }


      //---------------------[ assignment operations ]---------------------------------

      /**
      * Copy operator
      * @param rhs Right hand side argument of binary operator.
      */
      SIMD_INLINE IMatrix3x3<T>& operator=(const IMatrix3x3<T>& rhs)
      {
          std::memcpy(mData, rhs.mData, sizeof(T) * components);
          return *this;
      }

      /**
      * Copy casting operator
      * @param rhs Right hand side argument of binary operator.
      */
      template<class FromT>
      SIMD_INLINE IMatrix3x3<T>& operator=(const IMatrix3x3<FromT>& rhs)
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
      SIMD_INLINE IMatrix3x3<T>& operator=(const T* rhs)
      {
          std::memcpy(mData, rhs, sizeof(T) * components);
          return *this;
      }


      /// Overloaded operator for equality condition
      SIMD_INLINE bool operator == (const IMatrix3x3<T>& matrix) const
      {
          return (mRows[0] == matrix.mRows[0] &&
                  mRows[1] == matrix.mRows[1] &&
                  mRows[2] == matrix.mRows[2] );
      }

      /// Overloaded operator for the is different condition
      SIMD_INLINE bool operator != (const IMatrix3x3<T>& matrix) const
      {
            return !(*this == matrix);
      }




      //-------------[ conversion data ]-----------------------------

      /**
       * Conversion to pointer operator
       * @return Pointer to internally stored (in management of class IMatrix3x3<T>)
       * used for passing IMatrix3x3<T> values to gl*[fd]v functions.
       */
      SIMD_INLINE operator T*()
      {
          return &mRows[0][0];
      }

      /**
       * Conversion to pointer operator
       * @return Constant Pointer to internally stored (in management of class IMatrix3x3<T>)
       * used for passing IMatrix3x3<T> values to gl*[fd]v functions.
       */
      SIMD_INLINE operator const T*() const
      {
          return &mRows[0][0];
      }


      /**
      * Conversion to pointer operator
      */
      SIMD_INLINE const T* GetData() const
      {
          return mData;
      }

     /**
      * Conversion to pointer operator
      */
      SIMD_INLINE T* GetData()
      {
          return mData;
      }



      /// Overloaded operator to read element of the matrix.
      SIMD_INLINE const IVector3D<T>& operator[](int row) const
      {
           return mRows[row];
      }

      /// Overloaded operator to read/write element of the matrix.
      SIMD_INLINE IVector3D<T>& operator[](int row)
      {
           return mRows[row];
      }




      /// Return a column
      SIMD_INLINE IVector3D<T> GetColumn(int i) const
      {
          assert(i>= 0 && i<3);
          return IVector3D<T> (mRows[0][i], mRows[1][i], mRows[2][i]);
      }

      /// Return a row
      SIMD_INLINE IVector3D<T> GetRow(int i) const
      {
          assert(i>= 0 && i<3);
          return mRows[i];
      }




      //--------------------[ matrix with scalar operations ]---------------------
      /**
      * Addition operator
      * @param rhs Right hand side argument of binary operator.
      */
      SIMD_INLINE IMatrix3x3<T> operator+(T rhs) const
      {
          IMatrix3x3<T> ret;
          for (int i = 0; i < 3; i++) ret.mRows[i] = mRows[i] + rhs;
          return ret;
      }

      /**
      * Subtraction operator
      * @param rhs Right hand side argument of binary operator.
      */
      SIMD_INLINE IMatrix3x3<T> operator-(T rhs) const
      {
          IMatrix3x3<T> ret;
          for (int i = 0; i < 3; i++) ret.mRows[i] = mRows[i] - rhs;
          return ret;
      }


      //--------------------[ matrix with matrix operations ]---------------------



      /**
       * Addition operator
       * @param rhs Right hand side argument of binary operator.
       */
      SIMD_INLINE IMatrix3x3<T> operator+(const IMatrix3x3<T>& rhs) const
      {
          IMatrix3x3<T> ret;
    	  for (int i = 0; i < 3; i++)   ret.mRows[i] = mRows[i] + rhs.mRows[i];
    	  return ret;
      }

      /**
       * Subtraction operator
       * @param rhs Right hand side argument of binary operator.
       */
      SIMD_INLINE IMatrix3x3<T> operator-(const IMatrix3x3<T>& rhs) const
      {
          IMatrix3x3<T> ret;
    	  for (int i = 0; i < 3; i++)   ret.mRows[i] = mRows[i] - rhs.mRows[i];
    	  return ret;
      }

      //---------------------------- Friend method -------------------------//

      /// Overloaded operator for the negative of the matrix
      friend SIMD_INLINE IMatrix3x3<T>  operator-(const IMatrix3x3<T> & matrix)
      {
          return IMatrix3x3<T>(-matrix.mRows[0][0], -matrix.mRows[0][1], -matrix.mRows[0][2],
                               -matrix.mRows[1][0], -matrix.mRows[1][1], -matrix.mRows[1][2],
                               -matrix.mRows[2][0], -matrix.mRows[2][1], -matrix.mRows[2][2]);
      }

      /// Overloaded operator for multiplication with a number
      friend SIMD_INLINE IMatrix3x3<T>  operator*(T nb, const IMatrix3x3<T>& matrix)
      {
          return IMatrix3x3<T>(matrix.mRows[0][0] * nb, matrix.mRows[0][1] * nb, matrix.mRows[0][2] * nb,
                               matrix.mRows[1][0] * nb, matrix.mRows[1][1] * nb, matrix.mRows[1][2] * nb,
                               matrix.mRows[2][0] * nb, matrix.mRows[2][1] * nb, matrix.mRows[2][2] * nb);
      }

      /// Overloaded operator for multiplication with a matrix
      friend SIMD_INLINE IMatrix3x3<T>  operator*(const IMatrix3x3<T>& matrix, T nb)
      {
          return nb * matrix;
      }

      /// Overloaded operator for inveret multiplication with a number
      friend SIMD_INLINE IMatrix3x3<T>  operator/(T nb, const IMatrix3x3<T>& matrix)
      {
          return IMatrix3x3<T>(matrix.mRows[0][0] / nb, matrix.mRows[0][1] / nb, matrix.mRows[0][2] / nb,
                               matrix.mRows[1][0] / nb, matrix.mRows[1][1] / nb, matrix.mRows[1][2] / nb,
                               matrix.mRows[2][0] / nb, matrix.mRows[2][1] / nb, matrix.mRows[2][2] / nb);
      }

      /// Overloaded operator for inveret multiplication with a matrix
      friend SIMD_INLINE IMatrix3x3<T>  operator/(const IMatrix3x3<T>& matrix, T nb)
      {
          return nb / matrix;
      }


      //-------------------------------------------------------------------------//


      //-----------------------------------------------------------//

          SIMD_INLINE void Scale(T x, T y, T z)
          {
              *this = *this * CreateScale(x,y,z);
          }

          SIMD_INLINE void Scale(T factor)
          {
              *this = *this * CreateScale(factor);
          }

          SIMD_INLINE void Scale(const IVector3D<T>& scale)
          {
              *this = *this * CreateScale(scale);
          }

          //----------------------------------------------//

          SIMD_INLINE void Translate(T x, T y)
          {
              *this = *this * CreateTranslation(IVector2D<T>(x,y));
          }

          SIMD_INLINE void Translate(const IVector2D<T>& pos)
          {
              *this = *this * CreateTranslation(pos.x,pos.y);
          }

          //----------------------------------------------//

          SIMD_INLINE void Rotate( T x, T y, T z )
          {
              *this = *this * CreateEulerAnglesToRotationMatrix(IVector3D<T>(x,y,z));
          }

          SIMD_INLINE void RotateAxis(T angle, const IVector3D<T>& axis)
          {
              *this = *this * CreateRotationAxis(axis,angle);
          }

          SIMD_INLINE void Rotate(const IQuaternion<T>& quaternion)
          {
              *this = *this * CreateRotation(quaternion);
          }

          //-----------------------------------------------------------//


      //--------------------[ multiply operators ]--------------------------------


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
      friend SIMD_INLINE IVector3D<T> operator*(const IVector3D<T>& rhs , const IMatrix3x3<T>& lhs)
      {
          return IVector3D<T>(lhs.mRows[0][0]*rhs.x + lhs.mRows[1][0]*rhs.y + lhs.mRows[2][0]*rhs.z,
                              lhs.mRows[0][1]*rhs.x + lhs.mRows[1][1]*rhs.y + lhs.mRows[2][1]*rhs.z,
                              lhs.mRows[0][2]*rhs.x + lhs.mRows[1][2]*rhs.y + lhs.mRows[2][2]*rhs.z);
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
      SIMD_INLINE IVector3D<T> operator*(const IVector3D<T>& rhs) const
      {
          return IVector3D<T>(mRows[0][0]*rhs.x + mRows[0][1]*rhs.y + mRows[0][2]*rhs.z,
                              mRows[1][0]*rhs.x + mRows[1][1]*rhs.y + mRows[1][2]*rhs.z,
                              mRows[2][0]*rhs.x + mRows[2][1]*rhs.y + mRows[2][2]*rhs.z);

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
      SIMD_INLINE IMatrix3x3<T> operator*(IMatrix3x3<T> rhs) const
      {
          IMatrix3x3<T> w;
    	  for(int i = 0; i < 3; i++)
    	  {
    		  for (int j = 0; j < 3; j++)
    		  {
    			  T n = 0;
    			  for (int k = 0; k < 3; k++)
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
      * Return the determinant minor of the matrix
      */
      SIMD_INLINE T DeterminantOfMinor( int  theRowHeightY , int  theColumnWidthX ) const
      {
          int x1 = theColumnWidthX == 0 ? 1 : 0;  /* always either 0 or 1 */
          int x2 = theColumnWidthX == 2 ? 1 : 2;  /* always either 1 or 2 */
          int y1 = theRowHeightY   == 0 ? 1 : 0;  /* always either 0 or 1 */
          int y2 = theRowHeightY   == 2 ? 1 : 2;  /* always either 1 or 2 */

          return ( mRows[y1][x1]  *  mRows[y2][x2] )
              -  ( mRows[y1][x2]  *  mRows[y2][x1] );
      }

      /**
       * Computes determinant of matrix
       * @return Determinant of matrix
       * @note This function does.
       */
      SIMD_INLINE T Determinant() const
      {
          return ( mRows[0][0] * DeterminantOfMinor(0,0) )
              -  ( mRows[1][0] * DeterminantOfMinor(1,0) )
              +  ( mRows[2][0] * DeterminantOfMinor(2,0) );
      }

      /**
      * Computes inverse matrix
      * @return Inverse matrix of this matrix.
      * @note This is a little bit time consuming operation
      */
      SIMD_INLINE IMatrix3x3<T> Inverse() const
      {
          // Compute the determinant of the matrix
            T determinant = Determinant();

            // Check if the determinant is equal to zero
             assert(IAbs(determinant) > MACHINE_EPSILON);

            T invDeterminant = T(1.0) / determinant;

            IMatrix3x3<T> tempMatrix((mRows[1][1]*mRows[2][2]-mRows[2][1]*mRows[1][2]),
                                     -(mRows[0][1]*mRows[2][2]-mRows[2][1]*mRows[0][2]),
                                      (mRows[0][1]*mRows[1][2]-mRows[0][2]*mRows[1][1]),
                                     -(mRows[1][0]*mRows[2][2]-mRows[2][0]*mRows[1][2]),
                                      (mRows[0][0]*mRows[2][2]-mRows[2][0]*mRows[0][2]),
                                     -(mRows[0][0]*mRows[1][2]-mRows[1][0]*mRows[0][2]),
                                      (mRows[1][0]*mRows[2][1]-mRows[2][0]*mRows[1][1]),
                                     -(mRows[0][0]*mRows[2][1]-mRows[2][0]*mRows[0][1]),
                                      (mRows[0][0]*mRows[1][1]-mRows[0][1]*mRows[1][0]));

            // Return the inverse matrix
          return (invDeterminant * tempMatrix);
      }



      /**
      * Transpose matrix.
      */
      SIMD_INLINE IMatrix3x3<T> Transpose() const
      {
          IMatrix3x3<T> ret;
          for (int i = 0; i < 3; i++)
          {
              for (int j = 0; j < 3; j++)
              {
                  ret.mRows[i][j] = mRows[j][i];
              }
          }
          return ret;

      }


      /**
      * Return the matrix with absolute values
      */
      SIMD_INLINE IMatrix3x3<T> AbsoluteMatrix() const
      {
          return IMatrix3x3<T>(IAbs(mRows[0][0]), IAbs(mRows[0][1]), IAbs(mRows[0][2]),
                               IAbs(mRows[1][0]), IAbs(mRows[1][1]), IAbs(mRows[1][2]),
                               IAbs(mRows[2][0]), IAbs(mRows[2][1]), IAbs(mRows[2][2]));
      }


      /**
      * Return the trace of the matrix
      */
      SIMD_INLINE T Trace() const
      {
          // Compute and return the trace
          return (mRows[0][0] + mRows[1][1] + mRows[2][2]);
      }



      /**
      * Return the Euler angle of the matrix
      */
      SIMD_INLINE IVector3D<T> GetEulerAngles() const
      {
          T rotXangle = IAtan2(-GetRow(1).z,
                                GetRow(2).z);

          T cosYangle = ISqrt(pow(GetRow(0).x, 2) +
                              pow(GetRow(0).y, 2));

          T rotYangle = IAtan2(GetRow(0).z, cosYangle);

          T sinXangle = ISin(rotXangle);
          T cosXangle = ICos(rotXangle);
          T rotZangle = IAtan2(cosXangle * GetRow(1).x +
                               sinXangle * GetRow(2).x,

          cosXangle * GetRow(1).y +
          sinXangle * GetRow(2).y);

          return IVector3D<T>(rotXangle,
                              rotYangle,
                              rotZangle);
      }


//      SIMD_INLINE bool closeEnough(const T& a, const T& b, const T& epsilon = std::numeric_limits<T>::epsilon()) const
//      {
//          return (epsilon > std::abs(a - b));
//      }


//      /**
//      * Return the Euler angle of the matrix
//      */
//      SIMD_INLINE IVector3D<T> GetEulerAngles() const
//      {
//          const IMatrix3x3<T> R(*this);
//          //check for gimbal lock
//          if (closeEnough(R[0][2], T(-1.0f)))
//          {
//              T x = 0; //gimbal lock, value of x doesn't matter
//              T y = M_PI / 2;
//              T z = x + atan2(R[1][0], R[2][0]);
//              return { x, y, z };
//          }
//          else if (closeEnough(R[0][2], 1.0f))
//          {
//              T x = 0;
//              T y = -M_PI / 2;
//              T z = -x + atan2(-R[1][0], -R[2][0]);
//              return { x, y, z };
//          }
//          else
//          {
//              //two solutions exist
//              T x1 = -asin(R[0][2]);
//              T x2 = M_PI - x1;

//              T y1 = atan2(R[1][2] / cos(x1), R[2][2] / cos(x1));
//              T y2 = atan2(R[1][2] / cos(x2), R[2][2] / cos(x2));

//              T z1 = atan2(R[0][1] / cos(x1), R[0][0] / cos(x1));
//              T z2 = atan2(R[0][1] / cos(x2), R[0][0] / cos(x2));

//              //choose one solution to return
//              //for example the "shortest" rotation
//              if ((std::abs(x1) + std::abs(y1) + std::abs(z1)) <=
//                  (std::abs(x2) + std::abs(y2) + std::abs(z2)))
//              {
//                  return { x1, y1, z1 };
//              }
//              else
//              {
//                  return { x2, y2, z2 };
//              }
//          }
//      }


      /**
      * Return the diagonalize of the matrix
      */
      SIMD_INLINE IMatrix3x3<T> Diagonalize( T threshold, int maxSteps )
      {
          IMatrix3x3<T> rot;
          rot.SetToIdentity();
          for (int step = maxSteps; step > 0; step--)
          {
            // find off-diagonal element [p][q] with largest magnitude
            int p = 0;
            int q = 1;
            int r = 2;
            T max = IAbs(mRows[0][1]);
            T v = IAbs(mRows[0][2]);
            if (v > max)
            {
                q = 2;
                r = 1;
                max = v;
            }
            v = IAbs(mRows[1][2]);
            if (v > max)
            {
                p = 1;
                q = 2;
                r = 0;
                max = v;
            }

            T t = threshold * (IAbs(mRows[0][0]) + IAbs(mRows[1][1]) + IAbs(mRows[2][2]));
            if (max <= t)
            {
                if (max <= MACHINE_EPSILON * t)
                {
                    break;
                }
                step = 1;
            }

            // compute Jacobi rotation J which leads to a zero for element [p][q]
            T mpq = mRows[p][q];
            T theta = (mRows[q][q] - mRows[p][p]) / (2 * mpq);
            T theta2 = theta * theta;
            T cos;
            T sin;
            if (theta2 * theta2 < T(10 / MACHINE_EPSILON))
            {
                t = (theta >= 0) ? 1 / (theta + ISqrt(1 + theta2))
                                 : 1 / (theta - ISqrt(1 + theta2));
                cos = 1 / ISqrt(1 + t * t);
                sin = cos * t;
            }
            else
            {
                // approximation for large theta-value, i.e., a nearly diagonal matrix
                t = 1 / (theta * (2 + T(0.5) / theta2));
                cos = 1 - T(0.5) * t * t;
                sin = cos * t;
            }

            // apply rotation to matrix (this = J^T * this * J)
            mRows[p][q] = mRows[q][p] = 0;
            mRows[p][p] -= t * mpq;
            mRows[q][q] += t * mpq;
            T mrp = mRows[r][p];
            T mrq = mRows[r][q];
            mRows[r][p] = mRows[p][r] = cos * mrp - sin * mrq;
            mRows[r][q] = mRows[q][r] = cos * mrq + sin * mrp;
            // apply rotation to rot (rot = rot * J)
            for (int i = 0; i < 3; i++)
            {
                IVector3D<T>& row = rot[i];
                mrp = row[p];
                mrq = row[q];
                row[p] = cos * mrp - sin * mrq;
                row[q] = cos * mrq + sin * mrp;
            }

          }

          return rot;
      }


      //--------------------------------------------------------------------//


      /// Overloaded operator for addition with assignment
      SIMD_INLINE IMatrix3x3<T>& operator+=(const IMatrix3x3<T>& matrix)
      {
          mRows[0] += matrix.mRows[0];
          mRows[1] += matrix.mRows[1];
          mRows[2] += matrix.mRows[2];
          return *this;
      }

      /// Overloaded operator for substraction with assignment
      SIMD_INLINE IMatrix3x3<T>& operator-=(const IMatrix3x3<T>& matrix)
      {
          mRows[0] -= matrix.mRows[0];
          mRows[1] -= matrix.mRows[1];
          mRows[2] -= matrix.mRows[2];
          return *this;
      }

      /// Overloaded operator for multiplication with a number with assignment
      SIMD_INLINE IMatrix3x3<T>& operator*=(T nb)
      {
          mRows[0] *= nb;
          mRows[1] *= nb;
          mRows[2] *= nb;
          return *this;
      }

      /// Overloaded operator for invert multiplication with a number with assignment
      SIMD_INLINE IMatrix3x3<T> &operator/=(T nb)
      {
          mRows[0] /= nb;
          mRows[1] /= nb;
          mRows[2] /= nb;
          return *this;
      }





      ///---------------------------[ Pulgins ] -----------------------------------///


      /**
       * @brief GrammSchmidt
       * @param dir
       * @return
       * Gramm schmidt process
       * https://www.math.hmc.edu/calculus/tutorials/gramschmidt/gramschmidt.pdf
       */
      static SIMD_INLINE IMatrix3x3<T> GrammSchmidt(const IVector3D<T>& dir)
      {
          IVector3D<T> up = IVector3D<T>::ZERO;
          IVector3D<T> right = IVector3D<T>::ZERO;
          IVector3D<T> front (dir);

          front = front * (1.0f / ISqrt(front.Dot(front)));
          if (IAbs(front.z) > T(0.5 + 0.001))
          {
              right = front.Cross(IVector3D<T> (-front.y, front.z, T(0.0)));
          }
          else
          {
              right = front.Cross(IVector3D<T> (-front.y, front.x, T(0.0)));
          }

          right = right * (T(1.0f) / ISqrt(right.Dot(right)));
          up = right.Cross(front);

          IMatrix3x3<T> m;
          m.mRows[0] = front;
          m.mRows[1] = up;
          m.mRows[2] = right;

          return m;
      }

      /// Return a skew-symmetric matrix using a given vector that can be used
      /// to compute Cross product with another vector using matrix multiplication
      static SIMD_INLINE IMatrix3x3<T> ComputeSkewSymmetricMatrixForCrossProduct(const IVector3D<T>& vector)
      {
          return IMatrix3x3<T>(0, -vector.z, vector.y,
                               vector.z, 0, -vector.x,
                              -vector.y, vector.x, 0);
      }




       /**
        * Returns a scaling matrix that scales by `factor` uniformly.
        *
        * @param scale Uniform scale factor.
        * @return Scaling matrix.
        */
       static SIMD_INLINE IMatrix3x3<T> CreateScale(const T _scale)
       {
           static IMatrix3x3 res;

           T _x = T(0. + _scale);
           T _y = T(0. + _scale);
           T _z = T(0. + _scale);

           res.mRows[0] = IVector3D<T>(_x , 0.f, 0.f);
           res.mRows[1] = IVector3D<T>(0.f, _y , 0.f);
           res.mRows[2] = IVector3D<T>(0.f, 0.f, _z );

           return res;
       }


       /**
        * Returns a scaling matrix that scales by `scaleFactors.x` and
        * 'scaleFactors.y' in the x and y axes respectively.
        *
        * @param scaleFactors Scale factors.
        * @return Scaling matrix.
        */
       static SIMD_INLINE IMatrix3x3<T> CreateScale( T _x , T _y , T _z  )
       {
           static IMatrix3x3 res;

//           res.mRows[0] = IVector3D<T>(_x , 0.f, 0.f);
//           res.mRows[1] = IVector3D<T>(0.f, _y , 0.f);
//           res.mRows[2] = IVector3D<T>(0.f, 0.f, _z );

           res.mRows[0] = IVector3D<T>(0.f, 0.f,  _z);
           res.mRows[1] = IVector3D<T>(0.f, _y , 0.f);
           res.mRows[2] = IVector3D<T>(_x , 0.f, 0.f);



           return res;
       }


       /**
        * Returns a scaling matrix that scales by `scaleFactors.x` and
        * 'scaleFactors.y' in the x and y axes respectively.
        *
        * @param scaleFactors Scale factors.
        * @return Scaling matrix.
        */
       static SIMD_INLINE IMatrix3x3<T> CreateScale(const IVector3D<T>& _scale )
       {
           static IMatrix3x3 res;

           T _x = T(0. + _scale.x);
           T _y = T(0. + _scale.y);
           T _z = T(0. + _scale.z);

           res.mRows[0] = IVector3D<T>(_x , 0.f, 0.f);
           res.mRows[1] = IVector3D<T>(0.f, _y , 0.f);
           res.mRows[2] = IVector3D<T>(0.f, 0.f, _z );

           return res;
       }



       /**
        * Returns a scaling around axis matrix that scales
        * @return axis to scaling matrix.
        */
       static SIMD_INLINE IMatrix3x3<T>  CreateScaleAroundAxis( const IVector3D<T> _axis , const T &_scale )
       {
           static IMatrix3x3<T> M;
           T bgamma = (_scale - 1.0);

           M[0][0] = 1.0+((bgamma)*((_axis.x * _axis.x)));
           M[1][0] =     ((bgamma)*((_axis.y * _axis.x)));
           M[2][0] =     ((bgamma)*((_axis.z * _axis.x)));

           M[0][1] =     ((bgamma)*((_axis.x * _axis.y)));
           M[1][1] = 1.0+((bgamma)*((_axis.y * _axis.y)));
           M[2][1] =     ((bgamma)*((_axis.z * _axis.y)));

           M[0][2] =      ((bgamma)*((_axis.x * _axis.z)));
           M[1][2] =      ((bgamma)*((_axis.y * _axis.z)));
           M[2][2] =  1.0+((bgamma)*((_axis.z * _axis.z)));

           return M;
       }

      /*****************************************************
       *  Help info to web site:  https://arxiv.org/pdf/1103.0156.pdf
       *****************************************************/
       /// Return lorentz demission distance world
      static SIMD_INLINE IMatrix3x3<T> CreateLorentzRotationBoost(  const IVector3D<T>& vel , const T &_LightSpeed = DEFAUL_LIGHT_MAX_VELOCITY_C )
      {
          const IVector3D<T> n = vel.GetUnit();
          const T            v = vel.Length();

          static IMatrix3x3<T> M;

          const T c = _LightSpeed;

          //boost this Lorentz vector
          T gamma = 1.0 * ISqrt( 1.0 - (v*v) / (c*c) );

         // T bgamma = gamma * gamma / (1.0 + gamma);
          T bgamma = (gamma - 1.0);


          M[0][0] = 1.0+((bgamma)*((n.x * n.x)));
          M[1][0] =     ((bgamma)*((n.y * n.x)));
          M[2][0] =     ((bgamma)*((n.z * n.x)));


          M[0][1] =     ((bgamma)*((n.x * n.y)));
          M[1][1] = 1.0+((bgamma)*((n.y * n.y)));
          M[2][1] =     ((bgamma)*((n.z * n.y)));


          M[0][2] =     ((bgamma)*((n.x * n.z)));
          M[1][2] =     ((bgamma)*((n.y * n.z)));
          M[2][2] = 1.0+((bgamma)*((n.z * n.z)));

          return M;
      }


      static SIMD_INLINE IMatrix3x3<T> CreateLorentzRotationBoost(  const T &gamma , const IVector3D<T> n  )
      {

              static IMatrix3x3<T> M;

             // T bgamma = gamma * gamma / (1.0 + gamma);
              T bgamma = (gamma - 1.0);


              M[0][0] = 1.0+((bgamma)*((n.x * n.x)));
              M[1][0] =     ((bgamma)*((n.y * n.x)));
              M[2][0] =     ((bgamma)*((n.z * n.x)));


              M[0][1] =     ((bgamma)*((n.x * n.y)));
              M[1][1] = 1.0+((bgamma)*((n.y * n.y)));
              M[2][1] =     ((bgamma)*((n.z * n.y)));


              M[0][2] =     ((bgamma)*((n.x * n.z)));
              M[1][2] =     ((bgamma)*((n.y * n.z)));
              M[2][2] = 1.0+((bgamma)*((n.z * n.z)));

              return M;
       }


      /// Creates translation matrix
      /**
       * Creates translation matrix.
       * @param x X-direction translation
       * @param y Y-direction translation
       * @param z for Z-coordinate translation (implicitly Set to 1)
       */
      static SIMD_INLINE IMatrix3x3<T> CreateTranslation(T x, T y, T z = 1)
      {
          IMatrix3x3<T> ret;
          ret.mRows[2][0] = x;
          ret.mRows[2][1] = y;
          ret.mRows[2][2] = z;

          return ret;
      }


      /// Creates translation matrix
      /**
       * Creates translation matrix.
       * @param x X-direction translation
       * @param y Y-direction translation
       * @param z for Z-coordinate translation (implicitly Set to 1)
       */
      static SIMD_INLINE IMatrix3x3<T> CreateTranslation(const IVector2D<T> & v, T z = 1)
      {
          IMatrix3x3<T> ret;
          ret.mRows[2][0] = v.x;
          ret.mRows[2][1] = v.y;
          ret.mRows[2][2] = z;

          return ret;
      }


      static SIMD_INLINE IMatrix3x3<T> CreateRotationEulerAngle(const T& Roll , const T& Pitch , const T& Yaw )
      {
          IMatrix3x3<T> M;

          const T SR = ISin(Roll  /** M_PI / 180.f*/);
          const T SP = ISin(Pitch /** M_PI / 180.f*/);
          const T SY = ISin(Yaw   /** M_PI / 180.f*/);
          const T CR = ICos(Roll  /** M_PI / 180.f*/);
          const T CP = ICos(Pitch /** M_PI / 180.f*/);
          const T CY = ICos(Yaw   /** M_PI / 180.f*/);

          M[0][0]	= CP * CY;
          M[0][1]	= CP * SY;
          M[0][2]	= SP;

          M[1][0]	= SR * SP * CY - CR * SY;
          M[1][1]	= SR * SP * SY + CR * CY;
          M[1][2]	= - SR * CP;

          M[2][0]	= -( CR * SP * CY + SR * SY );
          M[2][1]	= CY * SR - CR * SP * SY;
          M[2][2]	= CR * CP;

          return M;
      }


      // Return a 4x4 rotation of axis to matrix
      static SIMD_INLINE IMatrix3x3<T> CreateRotationAxis(const IVector3D<T>& axis, T angle)
      {

          //angle = angle / 180.0f * (T)M_PI;

          T cosA = ICos(angle);
          T sinA = ISin(angle);
          IMatrix3x3<T> rotationMatrix;

          rotationMatrix.mRows[0][0] = cosA + (1-cosA) * axis.x * axis.x;
          rotationMatrix.mRows[1][0] = (1-cosA) * axis.x * axis.y - axis.z * sinA;
          rotationMatrix.mRows[2][0] = (1-cosA) * axis.x * axis.z + axis.y * sinA;

          rotationMatrix.mRows[0][1] = (1-cosA) * axis.x * axis.y + axis.z * sinA;
          rotationMatrix.mRows[1][1] = cosA + (1-cosA) * axis.y * axis.y;
          rotationMatrix.mRows[2][1] = (1-cosA) * axis.y * axis.z - axis.x * sinA;


          rotationMatrix.mRows[0][2] = (1-cosA) * axis.x * axis.z - axis.y * sinA;
          rotationMatrix.mRows[1][2] = (1-cosA) * axis.y * axis.z + axis.x * sinA;
          rotationMatrix.mRows[2][2] = cosA + (1-cosA) * axis.z * axis.z;


          return rotationMatrix;
      }



      static SIMD_INLINE IMatrix3x3<T> CreateRotation(const IQuaternion<T>& Quat)
      {
          static IMatrix3x3<T> M;

          T D1, D2, D3, D4, D5, D6, D7, D8, D9; //Dummy variables to hold precalcs

          D1 = (Quat.v.x * Quat.v.x) * 2.0f;
          D2 = (Quat.v.y * Quat.v.y) * 2.0f;
          D3 = (Quat.v.z * Quat.v.z) * 2.0f;

          T RTimesTwo = Quat.w * 2.0f;
          D4 = Quat.v.x * RTimesTwo;
          D5 = Quat.v.y * RTimesTwo;
          D6 = Quat.v.z * RTimesTwo;

          D7 = (Quat.v.x * Quat.v.y) * 2.0f;
          D8 = (Quat.v.x * Quat.v.z) * 2.0f;
          D9 = (Quat.v.y * Quat.v.z) * 2.0f;

          M.mRows[0][0] = 1.0f - D2 - D3;
          M.mRows[1][0] = D7 - D6;
          M.mRows[2][0] = D8 + D5;

          M.mRows[0][1] = D7 + D6;
          M.mRows[1][1] = 1.0f - D1 - D3;
          M.mRows[2][1] = D9 - D4;

          M.mRows[0][2] = D8 - D5;
          M.mRows[1][2] = D9 + D4;
          M.mRows[2][2] = 1.0f - D1 - D2;

          return M;
      }



      /// <summary>
      /// Creates a left-handed spherical billboard that rotates around a specified object position.
      /// </summary>
      /// <param name="objectPosition">The position of the object around which the billboard will rotate.</param>
      /// <param name="cameraPosition">The position of the camera.</param>
      /// <param name="cameraUpVector">The up vector of the camera.</param>
      /// <param name="cameraForwardVector">The forward vector of the camera.</param>
      /// <param name="result">When the method completes, contains the created billboard Matrix3x3.</param>
      static IVector3D<T> BillboardLH(const IVector3D<T> &objectPosition, const IVector3D<T> &cameraPosition, const IVector3D<T> &cameraUpVector, const IVector3D<T> &cameraForwardVector)
      {
          IVector3D<T> crossed;
          IVector3D<T> final;
          IVector3D<T> difference = cameraPosition - objectPosition;

          T lengthSq = difference.LengthSquare();
          if (IsZero(lengthSq))
          {
              difference = -cameraForwardVector;
          }
          else
          {
              difference *= (T)(1.0 / ISqrt(lengthSq));
          }

          crossed = Cross(cameraUpVector, difference );
          crossed.Normalize();
          final   = Cross(difference, crossed);

          IMatrix3x3 result;
          result[0][0] = crossed.x;
          result[0][1] = crossed.y;
          result[0][2] = crossed.z;
          result[1][0] = final.x;
          result[1][1] = final.y;
          result[1][2] = final.z;
          result[2][0] = difference.x;
          result[2][1] = difference.y;
          result[2][2] = difference.z;

          return result;
      }



      /// <summary>
      /// Creates a right-handed spherical billboard that rotates around a specified object position.
      /// </summary>
      /// <param name="objectPosition">The position of the object around which the billboard will rotate.</param>
      /// <param name="cameraPosition">The position of the camera.</param>
      /// <param name="cameraUpVector">The up vector of the camera.</param>
      /// <param name="cameraForwardVector">The forward vector of the camera.</param>
      /// <param name="result">When the method completes, contains the created billboard Matrix3x3.</param>
      static IVector3D<T> BillboardRH(const IVector3D<T> &objectPosition, const IVector3D<T> &cameraPosition, const IVector3D<T> &cameraUpVector, const IVector3D<T> &cameraForwardVector)
      {
          IVector3D<T> crossed;
          IVector3D<T> final;
          IVector3D<T> difference = objectPosition - cameraPosition;

          T lengthSq = difference.LengthSquare();
          if (IsZero(lengthSq))
          {
              difference = -cameraForwardVector;
          }
          else
          {
              difference *= (T)(1.0 / ISqrt(lengthSq));
          }

          crossed = Cross(cameraUpVector, difference);
          crossed.Normalize();
          final   = Cross(difference, crossed);

          IMatrix3x3 result;
          result[0][0] = crossed.x;
          result[0][1] = crossed.y;
          result[0][2] = crossed.z;
          result[1][0] = final.x;
          result[1][1] = final.y;
          result[1][2] = final.z;
          result[2][0] = difference.x;
          result[2][1] = difference.y;
          result[2][2] = difference.z;

          return result;
      }



      /// <summary>
      /// Creates a left-handed, look-at Matrix3x3.
      /// </summary>
      /// <param name="eye">The position of the viewer's eye.</param>
      /// <param name="target">The camera look-at target.</param>
      /// <param name="up">The camera's up vector.</param>
      /// <param name="result">When the method completes, contains the created look-at Matrix3x3.</param>
      static SIMD_INLINE IMatrix3x3<T> LookAtLH(const IVector3D<T> &eye, const IVector3D<T> &target, const IVector3D<T> &up)
      {
          IVector3D<T>  xaxis, yaxis, zaxis;
          zaxis = target - eye;
          if (IAbs(zaxis.x) < MACHINE_EPSILON &&
              IAbs(zaxis.y) < MACHINE_EPSILON &&
              IAbs(zaxis.z) < MACHINE_EPSILON)
          {
              return IMatrix3x3<T>::IDENTITY;
          }
          zaxis.Normalize();

          xaxis = Cross(up,zaxis);
          xaxis.Normalize();

          yaxis = Cross(zaxis, xaxis);

          return IMatrix3x3<T> (xaxis.x, yaxis.x, zaxis.x,
                                xaxis.y, yaxis.y, zaxis.y,
                                xaxis.z, yaxis.z, zaxis.z);



      }


      /// <summary>
      /// Creates a right-handed, look-at Matrix3x3.
      /// </summary>
      /// <param name="eye">The position of the viewer's eye.</param>
      /// <param name="target">The camera look-at target.</param>
      /// <param name="up">The camera's up vector.</param>
      /// <param name="result">When the method completes, contains the created look-at Matrix3x3.</param>
      static SIMD_INLINE IMatrix3x3<T> LookAtRH(const IVector3D<T> & eye, const IVector3D<T> & target, const IVector3D<T> & up )
      {
          IVector3D<T>  xaxis, yaxis, zaxis;

          zaxis = (eye - target);
          if (IAbs(zaxis.x) < MACHINE_EPSILON &&
              IAbs(zaxis.y) < MACHINE_EPSILON &&
              IAbs(zaxis.z) < MACHINE_EPSILON)
          {
              return IMatrix3x3<T>::IDENTITY;
          }
          zaxis.Normalize();

          xaxis = Cross(up,zaxis);
          xaxis.Normalize();

          yaxis = Cross(zaxis, xaxis);

          return IMatrix3x3<T> (xaxis.x, yaxis.x, zaxis.x,
                                xaxis.y, yaxis.y, zaxis.y,
                                xaxis.z, yaxis.z, zaxis.z);
      }




      //----------[ output operator ]----------------------------
      /**
      * Output to stream operator
      * @param lhs Left hand side argument of operator (commonly ostream instance).
      * @param rhs Right hand side argument of operator.
      * @return Left hand side argument - the ostream object passed to operator.
      */
      friend std::ostream& operator <<(std::ostream& lhs, const IMatrix3x3<T>& rhs)
      {
          for (int i = 0; i < 3; i++)
          {
              lhs << "|\t";
              for (int j = 0; j < 3; j++)
              {
                  lhs << rhs[i][j] << "\t";
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
    static const IMatrix3x3<T> IDENTITY;

    /**
     * The additive identity matrix.
     */
    static const IMatrix3x3<T> ZERO;
};



template<class T> const IMatrix3x3<T> IMatrix3x3<T>::IDENTITY = IMatrix3x3<T>(1.0, 0.0, 0.0,
																				 0.0, 1.0, 0.0,
																				 0.0, 0.0, 1.0);

template<class T> const IMatrix3x3<T> IMatrix3x3<T>::ZERO = IMatrix3x3<T>(0.0, 0.0, 0.0,
																			 0.0, 0.0, 0.0,
																			 0.0, 0.0, 0.0);


template<class T>  SIMD_INLINE IMatrix3x3<T> operator ^ (const IVector3D<T> lhs , const IVector3D<T> rhs )
{
    return IMatrix3x3<T>( lhs.x * rhs.x , lhs.x * rhs.y, lhs.x * rhs.z,
                          lhs.y * rhs.x , lhs.y * rhs.y, lhs.y * rhs.z,
                          lhs.z * rhs.x , lhs.z * rhs.y, lhs.z * rhs.z);
}


//--------------------------------------
// Typedef shortcuts for Matrix3x3
//-------------------------------------

using IMatrix3x3r    = IMatrix3x3<Real>;
using IMatrix3x3f    = IMatrix3x3<float>;
using IMatrix3x3d    = IMatrix3x3<double>;
using IMatrix3x3i    = IMatrix3x3<std::int32_t>;
using IMatrix3x3ui   = IMatrix3x3<std::uint32_t>;
using IMatrix3x3b    = IMatrix3x3<std::int8_t>;
using IMatrix3x3ub   = IMatrix3x3<std::uint8_t>;



} /* namespace */


#endif /* IMATRIX3X3_H_ */
