/********************************************************************************
*
* IMatrix.h
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

#ifndef IMATRIX_H
#define IMATRIX_H


#include <cmath>
#include <cstring>
#include <algorithm>
#include <cstdint>
#include <initializer_list>
#include <assert.h>
#include <vector>

#include "IReal.h"
#include "IFunc.h"

#include "IVector.h"


namespace IMath
{

using namespace std;


template <typename T, std::size_t Rows, std::size_t Cols>
class IMatrix;

//============================================//
namespace
{
    template <typename T, std::size_t Rows, std::size_t Cols>
    struct MatrixDefaultInitializer
    {
        static void Initialize(IMatrix<T, Rows, Cols>& matrix)
        {
            matrix.Reset();
        }
    };

    template <typename T, std::size_t N>
    struct MatrixDefaultInitializer<T, N, N>
    {
        static void Initialize(IMatrix<T, N, N>& matrix)
        {
            matrix.SetToIdentity();
        }
    };

#define GS_ROW_MAJOR_STORAGE true

#ifdef GS_ROW_MAJOR_STORAGE
#   define GS_FOREACH_ROW_COL(r, c)             \
        for (std::size_t r = 0; r < Rows; ++r)  \
        for (std::size_t c = 0; c < Cols; ++c)
#else
#   define GS_FOREACH_ROW_COL(r, c)             \
        for (std::size_t c = 0; c < Cols; ++c)  \
        for (std::size_t r = 0; r < Rows; ++r)
#endif

    //======================================================//
    //! Internal class for implementation details.
    template <template <typename, std::size_t, std::size_t>
              class M,  typename T,
              std::size_t Rows,
              std::size_t Cols>
    class MatrixHelper
    {

            MatrixHelper() = delete;

        public:


            //friend bool Gs::Inverse<T, Rows>(Matrix<T, Rows, Cols>&, const Matrix<T, Rows, Cols>&);
            static std::vector<T> MatrixToArray(const M<T, Rows, Cols>& mat)
            {
                std::vector<T> vec(Rows*Cols);
                for (std::size_t r = 0, i = 0; r < Rows; ++r)
                {
                    for (std::size_t c = 0; c < Cols; ++c)
                        vec[i++] = mat(r, c);
                }

                return vec;
            }

            static T OrderedDeterminant(const std::vector<T>& mat, std::size_t order)
            {
                if (order == 1) return mat[0];
                std::vector<T> minorMat((order - 1)*(order - 1), T());

                T det = T(0.0);
                for (std::size_t i = 0; i < order; ++i)
                {
                    GetMinorMatrix(mat, minorMat, 0, i, order);
                    if (i % 2 == 1)
                    {
                        det -= mat[i] * OrderedDeterminant(minorMat, order - 1);
                    }
                    else
                    {
                        det += mat[i] * OrderedDeterminant(minorMat, order - 1);
                    }
                }

                return det;
            }



        private:

            static void GetMinorMatrix( const std::vector<T>& mat, std::vector<T>& minorMat,
                                        std::size_t row, std::size_t column,
                                        std::size_t order)
            {
                for (std::size_t r = 1, i = 0; r < order; ++r)
                {
                    if (r != row)
                    {
                        for (std::size_t c = 0, j = 0; c < order; ++c)
                        {
                            if (c != column)
                            {
                                minorMat[i*(order - 1) + j] = mat[r*order + c];
                                ++j;
                            }
                        }
                        ++i;
                    }
                }
            }

    };
    //======================================================//

}




template <typename T, std::size_t Rows, std::size_t Cols>
class IMatrix
{

    public:

        static_assert(Rows*Cols > 0, "matrices must consist of at least 1x1 elements");

        /* ----- Static members ----- */

        //! Number of rows of this matrix type.
        static const std::size_t rows       = Rows;

        //! Number of columns of this matrix type.
        static const std::size_t columns    = Cols;

        //! Number of scalar elements of this matrix type.
        static const std::size_t elements   = Rows*Cols;

        /* ----- Typenames ----- */

        //! Specifies the typename of the scalar components.
        using ScalarType        = T;

        //! Typename of this matrix type.
        using ThisType          = IMatrix<T, Rows, Cols>;

        //! Typename of the transposed of this matrix type.
        using TransposedType    = IMatrix<T, Cols, Rows>;




     private:

         union
         {
            T m_[ThisType::elements];

#ifdef GS_ROW_MAJOR_STORAGE
            IVector<T,Rows>  mVec[Cols];
#else
            IVector<T,Cols>  mVec[Rows];
#endif
         };


     public:

         //! Deffault constructor.
         IMatrix()
         {
#ifndef GS_DISABLE_AUTO_INIT
             MatrixDefaultInitializer<T, Rows, Cols>::Initialize(*this);
#endif
         }


         //! Copy constructor.
         IMatrix(const ThisType& rhs)
         {
             *this = rhs;
         }


         //template <typename Type>
         IMatrix(const std::initializer_list<T>& values)
         {
             std::size_t i = 0, n = values.size();
             for (auto it = values.begin(); i < n; ++i, ++it)
             {
                 m_[i] = *it;
             }
         }


//         template <typename... Type>
//         IMatrix(Type&&... a)
//         {
//             IMatrix({std::forward<Type>(a)...});
//         }

         //! Initializes this matrix with the specified values (row by row, and column by column).
//         IMatrix(const std::initializer_list<T>& values)
//         {
//             std::size_t i = 0, n = values.size();
//             for (auto it = values.begin(); i < n; ++i, ++it)
//             {
//                 (*this)(i / columns, i % columns) = *it;
//             }

//             for (; i < elements; ++i)
//             {
//                 (*this)(i / columns, i % columns) = T(0);
//             }
//         }

         template <typename... Type>
         void Set(Type&&... a)
         {
             Set({std::forward<T>(a)...});
         }

         template <typename Type>
         void Set(std::initializer_list<Type>&& values)
         {
             std::size_t i = 0, n = values.size();
             for (auto it = values.begin(); i < n; ++i, ++it)
             {
                 m_[i] = (*it);
             }
         }


         //! Restes all matrix elements to zero.
         void SetToZero()
         {
             for (std::size_t i = 0; i < ThisType::elements; ++i)
             {
                 m_[i] = T(0);
             }
         }

         //! Loads the identity for this matrix.
         void SetToIdentity()
         {
             GS_FOREACH_ROW_COL(r, c)
             {
                 (*this)(r, c) = (r == c ? T(1.0) : T(0.0));
             }
         }

         //! Returns an identity matrix.
         static ThisType Identity()
         {
             ThisType result;
             result.LoadIdentity();
             return result;
         }


        /**
        \brief Returns a reference to a single matrix element at the specified location.
        \param[in] row Specifies the zero-based row index. Must be in the range [0, Rows).
        \param[in] col Specifies the zero-based column index. Must be in the range [0, Cols).
        \throws std::runtime_error If the macro 'GS_ENABLE_ASSERT' and the macro 'GS_ASSERT_EXCEPTION' are defined,
         and either the row or the column is out of range.
         **/
         T& operator () (std::size_t row, std::size_t col)
         {
             assert(row < Rows);
             assert(col < Cols);
#ifdef GS_ROW_MAJOR_STORAGE
             return m_[row*Cols + col];
#else
             return m_[col*Rows + row];
#endif
         }


         /**
         \brief Returns a constant reference to a single matrix element at the specified location.
         \param[in] row Specifies the zero-based row index. Must be in the range [0, Rows).
         \param[in] col Specifies the zero-based column index. Must be in the range [0, Cols).
         \throws std::runtime_error If the macro 'GS_ENABLE_ASSERT' and the macro 'GS_ASSERT_EXCEPTION' are defined,
          and either the row or the column is out of range.
          */
         const T& operator () (std::size_t row, std::size_t col) const
         {
             assert(row < Rows);
             assert(col < Cols);
#ifdef GS_ROW_MAJOR_STORAGE
             return m_[row*Cols + col];
#else
             return m_[col*Rows + row];
#endif
         }


#ifdef GS_ROW_MAJOR_STORAGE
         IVector<T,Cols>& operator [] (std::size_t index)
         {
             assert(index < ThisType::rows);
             return mVec[index];
         }

         const IVector<T,Cols>& operator [] (std::size_t index) const
         {
             assert(index < ThisType::rows);
             return mVec[index];
         }
#else
         IVector<T,Rows>& operator [] (std::size_t index)
         {
             assert(index < ThisType::cols);
             return mVec[index];
         }

         const IVector<T,Rows>& operator [] (std::size_t index) const
         {
             assert(index < ThisType::cols);
             return mVec[index];
         }
#endif

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
         friend SIMD_INLINE IVector<T,Cols> operator*(const IVector<T,Cols>& rhs , const ThisType& lhs)
         {
             IVector<T,Cols> v_res;
             for (std::size_t i = 0; i < IVector<T,Cols>::components; ++i)
             {
                 T res = 0;
                 for (std::size_t j = 0; j < IVector<T,Cols>::components; ++j)
                 {
                        res += lhs.At(j,i) * rhs[j];
                 }
                 v_res.At(i) = res;
             }
             return v_res;
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
         SIMD_INLINE IVector<T,Cols> operator*(const IVector<T,Cols>& rhs) const
         {
             IVector<T,Cols> v_res;
             for (std::size_t i = 0; i < IVector<T,Cols>::components; ++i)
             {
                 T res = 0;
                 for (std::size_t j = 0; j < IVector<T,Cols>::components; ++j)
                 {
                        res += At(i,j) * rhs[j];
                 }
                 v_res.At(i) = res;
             }
             return v_res;
         }

         /**
         * Return the matrix with absolute values
         */
         SIMD_INLINE ThisType AbsoluteMatrix()// const
         {
             for (std::size_t i = 0; i < ThisType::elements; ++i)
                 m_[i] = IMath::IAbs(m_[i]);
             return *this;
         }

         SIMD_INLINE ThisType& operator += (const ThisType& rhs)
         {
             for (std::size_t i = 0; i < ThisType::elements; ++i)
                 m_[i] += rhs.m_[i];
             return *this;
         }

         SIMD_INLINE ThisType& operator -= (const ThisType& rhs)
         {
             for (std::size_t i = 0; i < ThisType::elements; ++i)
                 m_[i] -= rhs.m_[i];
             return *this;
         }

         SIMD_INLINE ThisType& operator *= (const ThisType& rhs)
         {
             *this = (*this * rhs);
              return *this;
         }

         SIMD_INLINE ThisType& operator *= (const T& rhs)
         {
             for (std::size_t i = 0; i < ThisType::elements; ++i)
             {
                 m_[i] *= rhs;
             }
             return *this;
         }

         SIMD_INLINE ThisType& operator /= (const T& rhs)
         {
             for (std::size_t i = 0; i < ThisType::elements; ++i)
             {
                 m_[i] /= rhs;
             }
             return *this;
         }

         SIMD_INLINE ThisType& operator = (const ThisType& rhs)
         {
             for (std::size_t i = 0; i < ThisType::elements; ++i)
             {
                 m_[i] = rhs.m_[i];
             }
             return *this;
         }

         //! Returns a transposed copy of this matrix.
         SIMD_INLINE TransposedType Transpose() const
         {
             TransposedType result;
             GS_FOREACH_ROW_COL(r, c)
             {
                 result(c, r) = (*this)(r, c);
             }
             return result;
         }

#ifdef GS_ROW_VECTORS

         T& At(std::size_t col, std::size_t row)
         {
             return (*this)(row, col);
         }

         const T& At(std::size_t col, std::size_t row) const
         {
             return (*this)(row, col);
         }




         /**
          \brief Returns the determinant of this matrix.
          \see Gs::Determinant
         */
         T Determinant() const
         {
             return Determinante(*this);
         }


         ThisType findMinor(const ThisType& pmatrix, std::size_t i, std::size_t j) const
         {
             ThisType minor;
             minor.LoadIdentity();

             std::size_t b = 0;
             std::size_t c = 0;
             for (std::size_t k = 0; k < Rows; k++)
             {
                 if (k == i) continue;
                 c = 0;
                 for (std::size_t s = 0; s < Cols; s++)
                 {
                     if (s == j) continue;
                     minor.At(b, c) = pmatrix.At(k, s);
                     c++;
                 }
                 b++;
             }

             return(minor);
         }


         ThisType Inverse() const
         {
             T k = Determinant();
             ThisType result = ThisType(*this);
             for (std::size_t i = 0; i < Cols; i++)
             {
                 for (std::size_t j = 0; j < Rows; j++)
                 {
                     ThisType minor(findMinor(*this, j, i));
                     result.At(i,j) = pow((-1), i + j) * minor.Determinant() / k;
                 }
             }
             return(result);
         }



         T Trace() const
         {
             static_assert(Rows == Cols, "traces can only be computed for squared matrices");

             T trace = T(0);
             for (std::size_t i = 0; i < Cols; ++i)
             {
                 trace += (*this)(i, i);
             }

             return trace;
         }


#else

         T& At(std::size_t row, std::size_t col)
         {
             return (*this)(row, col);
         }

         const T& At(std::size_t row, std::size_t col) const
         {
             return (*this)(row, col);
         }


         /**
          \brief Returns the determinant of this matrix.
          \see Gs::Determinant
         */
         SIMD_INLINE T Determinant() const
         {
             return Determinante(*this);
         }


         SIMD_INLINE ThisType findMinor(const ThisType& pmatrix, std::size_t i, std::size_t j) const
         {
             ThisType minor;
             minor.SetToIdentity();

             std::size_t b = 0;
             std::size_t c = 0;
             for (std::size_t k = 0; k < Rows; k++)
             {
                 if (k == i) continue;
                 c = 0;
                 for (std::size_t s = 0; s < Cols; s++)
                 {
                     if (s == j) continue;
                     minor.At(b, c) = pmatrix.At(k, s);
                     c++;
                 }
                 b++;
             }

             return(minor);
         }


         SIMD_INLINE ThisType Inverse() const
         {
             T k = Determinant();
             ThisType result = ThisType(*this);
             for (std::size_t i = 0; i < Rows; i++)
             {
                 for (std::size_t j = 0; j < Cols; j++)
                 {
                     ThisType minor(findMinor(*this, j, i));
                     result.At(i,j) = pow((-1), i + j) * minor.Determinant() / k;
                 }
             }
             return(result);
         }



         SIMD_INLINE T Trace() const
         {
             static_assert(Rows == Cols, "traces can only be computed for squared matrices");

             T trace = T(0);
             for (std::size_t i = 0; i < Rows; ++i)
             {
                 trace += (*this)(i, i);
             }

             return trace;
         }


         SIMD_INLINE ThisType Diagonalize()
         {
             ThisType diag_matrix(*this);	// Copy of our matrix (We must not crash our old matrix)
             std::size_t used_rows = 0;	// Iteratior for diag_matrix (the number of independent row)

             for (std::size_t j = 0; j < Cols; j++)
             {
                 bool row = false;
                 int i;
                 for (i = used_rows; i < Rows; i++) // (***)
                 {
                     if (diag_matrix.At(i,j) != 0)
                     {
                         row = true;
                         break;
                     }
                 }

                 if (!row)
                     continue;

                 T first_item = diag_matrix.At(i,j);

                 /*	We should select used rows.
                  Therefore we replace [i] and [used_rows] rows
                  and when we try to find new linery
                  independent row in (***) we do not consider
                  old rows by this technique.
                 */
                 for (std::size_t s = 0; s < Rows; s++)
                 {
                     T item_i = diag_matrix.At(i,s);
                     diag_matrix.At(i,s) = diag_matrix.At(used_rows,s);
                     diag_matrix.At(used_rows,s) = item_i;
                 }


                 used_rows++;

                 for (std::size_t s = used_rows; s < Rows; s++)
                 {
                     T koef = diag_matrix.At(s,j) / first_item;
                     for (int k = 0; k < Cols; k++)
                     {
                         diag_matrix.At(s,k) -= diag_matrix.At((used_rows - 1),k) * koef;
                     }
                 }
             }

             for (std::size_t i = used_rows+1; i < Rows; i++) // All other rows is linearly dependent.
             {
                 for (std::size_t j = 0; j < Cols; j++)
                 {
                     diag_matrix.At(i,j) = 0;
                 }
             }

             return diag_matrix;
         }

#endif


               //! Returns a pointer to the first element of this matrix.
               T* Ptr()
               {
                   return &(m_[0]);
               }

               //! Returns a constant pointer to the first element of this matrix.
               const T* Ptr() const
               {
                   return &(m_[0]);
               }

               /**
               Returns a type casted instance of this matrix.
               \tparam C Specifies the static cast type.
               */
               template <typename C> IMatrix<C, Rows, Cols> Cast() const
               {
                   IMatrix<C, Rows, Cols> result;

                   for (std::size_t i = 0; i < ThisType::elements; ++i)
                   {
                       result[i] = static_cast<C>(m_[i]);
                   }

                   return result;
               }


               /**
                * Equality operator.
                *
                * @param A The first matrix.
                * @param B The second matrix.
                * @return True if the two supplied matrices are equal. False otherwise.
                */
               friend SIMD_INLINE bool operator==(const ThisType& A, const ThisType& B)
               {
                   const T epsilon = MACHINE_EPSILON;
                   for (std::size_t i = 0; i < ThisType::elements; ++i)
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
               friend SIMD_INLINE  bool operator!=(const ThisType& A, const ThisType& B)
               {
                   return !(A == B);
               }


               //----------[ output operator ]----------------------------
               /**
               * Output to stream operator
               * @param lhs Left hand side argument of operator (commonly ostream instance).
               * @param rhs Right hand side argument of operator.
               * @return Left hand side argument - the ostream object passed to operator.
               */
               friend std::ostream& operator <<(std::ostream& lhs, const IMatrix<T,Rows,Cols>& rhs)
               {
                   for (int i = 0; i < Rows; i++)
                   {
                       lhs << "|\t";
                       for (int j = 0; j < Cols; j++)
                       {
                           lhs << rhs(i,j)  << "\t";
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

};


/* --- Global Operators --- */

template <typename T, std::size_t Rows, std::size_t Cols>
IMatrix<T, Rows, Cols> operator + (const IMatrix<T, Rows, Cols>& lhs, const IMatrix<T, Rows, Cols>& rhs)
{
    auto result = lhs;
    result += rhs;
    return result;
}

template <typename T, std::size_t Rows, std::size_t Cols>
IMatrix<T, Rows, Cols> operator - (const IMatrix<T, Rows, Cols>& lhs, const IMatrix<T, Rows, Cols>& rhs)
{
    auto result = lhs;
    result -= rhs;
    return result;
}

template <typename T, std::size_t Rows, std::size_t Cols>
IMatrix<T, Rows, Cols> operator * (const IMatrix<T, Rows, Cols>& lhs, const T& rhs)
{
    auto result = lhs;
    result *= rhs;
    return result;
}

template <typename T, std::size_t Rows, std::size_t Cols>
IMatrix<T, Rows, Cols> operator * (const T& lhs, const IMatrix<T, Rows, Cols>& rhs)
{
    auto result = rhs;
    result *= lhs;
    return result;
}

template <typename T, std::size_t Rows, std::size_t ColsRows, std::size_t Cols>
IMatrix<T, Rows, Cols> operator * (const IMatrix<T, Rows, ColsRows>& lhs,
                                   const IMatrix<T, ColsRows, Cols>& rhs)
{
    IMatrix<T, Rows, Cols> result;
    GS_FOREACH_ROW_COL(r, c)
    {
        result(c, r) = T(0.0);
        for (std::size_t i = 0; i < ColsRows; ++i)
        {
            result(c, r) += lhs(i, r)*rhs(c, i);
        }
    }
    return result;
}



/**
\brief Computes the determinant of an arbitrary NxN matrix.
\tparam M Specifies the matrix type. This should be "Matrix".
\tparam T Specifies the data type. This should be float or double.
\tparam Rows Specifies the rows of the matrix.
\tparam Cols Specifies the columns of the matrix.
\remarks The template arguments 'Rows' and 'Cols' must be equal, otherwise a compile time error will occur,
since a determinant is only defined for squared matrices.
\param[in] m Specifies the squared matrix for which the determinant is to be computed.
*/
template <typename T, std::size_t N>
T Determinante(const IMatrix<T, N, N>& m)
{
   using Helper = MatrixHelper<IMatrix, T, N, N>;
   return T(Helper::OrderedDeterminant(Helper::MatrixToArray(m), N));
}


//template <typename T, std::size_t N>
//bool Inverse(IMatrix<T, N, N>& inv, const IMatrix<T, N, N>& m)
//{
//    using Helper = MatrixHelper<IMatrix, T, N, N>;
//    return Helper::OrderedInverse(Helper::MatrixToArray(m) , N);
//    return false;//!!!
//}

//============================================//


} // /namespace Gs



#endif // IMATRIX_H
