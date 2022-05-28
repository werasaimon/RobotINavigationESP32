/********************************************************************************
 *
 * IQuaternion.h
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


#ifndef IQUATERNION_H_
#define IQUATERNION_H_


// Libraries
#include <cmath>
#include "IMatrix4x4.h"

namespace IMath
{

    /**
 * Quaternion class implementing some quaternion algebra operations.
 * Quaternion is kind of complex number it consists of its real part (w)
 * and its complex part v. This complex part has three elements, so we
 * can express it as xi + yj + zk . Note that coordinates of (x,y,z) are
 * hold inside v field.
 */
    template<class T> class IQuaternion
    {


     public:

        static_assert(std::is_floating_point<T>::value, "quaternions can only be used with floating point types");

        //! Specifies the typename of the scalar components.
        using ScalarType = T;

        //! Specifies the number of quaternion components. This is just for the internal template interface.
        static const std::size_t components = 4;


    private:


        // -------------------- Attributes -------------------- //

        /**
        * Real part of IQuaternion.
        **/
        T w;

        /**
        * Imaginary part of IQuaternion.
        **/
        union
        {
                struct
                {
                        T x;
                        T y;
                        T z;
                };

                IVector3D<T> v;
        };




        //---------------------- Help function -------------------//

        void twoaxisrot(T r11, T r12, T r21, T r31, T r32, T res[])
        {
          res[2] = IAtan2( r11, r12 );
          res[1] = IACos( r21 );
          res[0] = IAtan2( r31, r32 );
        }

        void threeaxisrot(T r11, T r12, T r21, T r31, T r32, T res[])
        {
          res[2] = IAtan2( r31, r32 );
          res[1] = IASin( r21 );
          res[0] = IAtan2( r11, r12 );
        }




    public:

        //--------------------------[ constructors ]------------------------------- //

        /**
        * IQuaternion constructor, Sets IQuaternion to (0 + 0i + 0j + 0k).
        */
        SIMD_INLINE IQuaternion()
        : w(0), x(0), y(0), z(0)
        {
        }

        /**
        * Copy constructor.
        */
        SIMD_INLINE IQuaternion(const IQuaternion<T>& q)
        : w(q.w), x(q.x) , y(q.y) , z(q.z)
        {
        }

        /**
        * Copy casting constructor.
        */
        template<class FromT>
        SIMD_INLINE IQuaternion(const IQuaternion<FromT>& q)
        : w(static_cast<T>(q.w)),
          x(static_cast<T>(q.x)),
          y(static_cast<T>(q.y)),
          z(static_cast<T>(q.z))
        {
        }

        /**
        * Creates IQuaternion object from real part w_ and complex part v_.
        * @param w_ Real part of IQuaternion.
        * @param v_ Complex part of IQuaternion (xi + yj + zk).
        */
        SIMD_INLINE IQuaternion(T w_, const IVector3D<T>& v)
        : w(w_), x(v.x) , y(v.y) , z(v.z)
        {
        }

        /**
        * Creates IQuaternion object from real part w_ and complex part v_.
        * @param w_ Real part of IQuaternion.
        * @param v_ Complex part of IQuaternion (xi + yj + zk).
        */
        SIMD_INLINE IQuaternion(const IVector3D<T>& v, T w_)
        : w(w_), x(v.x) , y(v.y) , z(v.z)
        {
        }

        /**
        * Creates IQuaternion object from value (w_ + xi + yj + zk).
        * @param x Complex coefficient for i complex constant.
        * @param y Complex coefficient for j complex constant.
        * @param z Complex coefficient for k complex constant.
        * @param w_ Real part of IQuaternion.
        */
        SIMD_INLINE IQuaternion(T _w,T _x, T _y, T _z)
        : w(_w), x(_x), y(_y), z(_z)
        {
        }


        /// Constructor which convert Euler angles (in radians) to a quaternion
        SIMD_INLINE IQuaternion(T angleX, T angleY, T angleZ)
        {
            SetEulerAngles(angleX, angleY, angleZ);
        }

        /// Constructor which convert Euler angles (in radians) to a quaternion
        SIMD_INLINE IQuaternion(const IVector3D<T>& eulerAngles)
        {
            SetEulerAngles(eulerAngles.x, eulerAngles.y, eulerAngles.z);
        }


        /// Constructor which convert Vector4D (in radians) to a quaternion
        SIMD_INLINE IQuaternion(const IVector4D<T>& _v4)
        {
            x = _v4.x;
            y = _v4.y;
            z = _v4.z;
            w = _v4.w;
        }

        /// Create a unit quaternion from a rotation matrix
        SIMD_INLINE IQuaternion(const IMatrix3x3<T>& matrix)
        {

            T trace = matrix[0].x + matrix[1].y + matrix[2].z;
            T temp[4];

            if (trace > T(0.0))
            {
                T s = ISqrt(trace + T(1.0));
                temp[3]=(s * T(0.5));
                s = T(0.5) / s;

                temp[0]=((matrix[2].y - matrix[1].z) * s);
                temp[1]=((matrix[0].z - matrix[2].x) * s);
                temp[2]=((matrix[1].x - matrix[0].y) * s);
            }
            else
            {
                int i = matrix[0].x < matrix[1].y ?
                (matrix[1].y < matrix[2].z ? 2 : 1) :
                (matrix[0].x < matrix[2].z ? 2 : 0);

                int j = (i + 1) % 3;
                int k = (i + 2) % 3;

                T s = ISqrt(matrix[i][i] - matrix[j][j] - matrix[k][k] + T(1.0));
                temp[i] = s * T(0.5);
                s = T(0.5) / s;

                temp[3] = (matrix[k][j] - matrix[j][k]) * s;
                temp[j] = (matrix[j][i] + matrix[i][j]) * s;
                temp[k] = (matrix[k][i] + matrix[i][k]) * s;
            }

            SetAllValues(temp[0],temp[1],temp[2],temp[3]);

        }

        /// Create a unit quaternion from a transform matrix
        SIMD_INLINE IQuaternion(const IMatrix4x4<T>& matrix)
        : IQuaternion( matrix.GetRotMatrix() )
        {

        }


        /**
        * Converts IQuaternion into rotation matrix.
        * @return Rotation matrix expressing this IQuaternion.
        */
        SIMD_INLINE IMatrix3x3<T> GetRotMatrix() const
        {
            IMatrix3x3<T> ret;

            T nQ = x*x + y*y + z*z + w*w;
            T s = 0.0;

            if(nQ > 0.0)
            {
                s = T(2.0) / nQ; /// normalize quaternion
            }

            // Computations used for optimization (less multiplications)
            T xs  = x*s;
            T ys  = y*s;
            T zs  = z*s;
            T wxs = w*xs;
            T wys = w*ys;
            T wzs = w*zs;
            T xxs = x*xs;
            T xys = x*ys;
            T xzs = x*zs;
            T yys = y*ys;
            T yzs = y*zs;
            T zzs = z*zs;

            // Create the matrix corresponding to the quaternion
            ret = IMatrix3x3<T>( T(1.0) - yys - zzs, xys-wzs, xzs + wys,
                                 xys + wzs, T(1.0) - xxs - zzs, yzs-wxs,
                                 xzs-wys, yzs + wxs, T(1.0) - xxs - yys );


            return ret;
        }



        /**
        * Converts IQuaternion into transformation matrix.
        * @note This method performs same operation as rotMatrix()
        * conversion method. But returns Matrix of 4x4 elements.
        * @return Transformation matrix expressing this IQuaternion.
        */
        SIMD_INLINE IMatrix4x4<T> TransformMatrix() const
        {
            IMatrix4x4<T> ret;

            T nQ = x*x + y*y + z*z + w*w;
            T s = 0.0;

            if(nQ > 0.0)
            {
                s = T(2.0) / nQ; /// normalize quaternion
            }

            // Computations used for optimization (less multiplications)
            T xs  = x*s;
            T ys  = y*s;
            T zs  = z*s;
            T wxs = w*xs;
            T wys = w*ys;
            T wzs = w*zs;
            T xxs = x*xs;
            T xys = x*ys;
            T xzs = x*zs;
            T yys = y*ys;
            T yzs = y*zs;
            T zzs = z*zs;

            // Create the matrix corresponding to the quaternion
            ret = IMatrix4x4<T>( T(1.0) - yys - zzs, xys-wzs, xzs + wys, T(0.0),
                                 xys + wzs, T(1.0) - xxs - zzs, yzs-wxs, T(0.0),
                                 xzs-wys, yzs + wxs, T(1.0) - xxs - yys, T(0.0),
                                 T(0.0),T(0.0),T(0.0),T(1.0));

            return ret;

        }

        //---------------------- Methods ---------------------//


        /// Set all the values
        SIMD_INLINE void SetAllValues(T newX, T newY, T newZ, T newW)
        {
            x = newX;
            y = newY;
            z = newZ;
            w = newW;
        }

        /// Set vector3 v
        SIMD_INLINE void SetV( const IVector3D<T>& _v)
        {
            x = _v.x;
            y = _v.y;
            z = _v.z;
        }

        /// Set component
        SIMD_INLINE void SetW(T _w)
        {
            w = _w;
        }
        SIMD_INLINE void SetX(T _x)
        {
            x = _x;
        }
        SIMD_INLINE void SetY(T _y)
        {
            y = _y;
        }
        SIMD_INLINE void SetZ(T _z)
        {
            z = _z;
        }


        /// Set the quaternion to zero
        SIMD_INLINE void SetToZero()
        {
            x = 0;
            y = 0;
            z = 0;
            w = 0;
        }

        /// Set to the identity quaternion
        SIMD_INLINE void SetToIdentity()
        {
            x = 0;
            y = 0;
            z = 0;
            w = 1;
        }

        SIMD_INLINE T R_component_1() const { return(w); }
        SIMD_INLINE T R_component_2() const { return(x); }
        SIMD_INLINE T R_component_3() const { return(y); }
        SIMD_INLINE T R_component_4() const { return(z); }


        SIMD_INLINE IVector3D<T> GetImage() const { return IVector3D<T>(x,y,z); }
        SIMD_INLINE            T GetReal() const { return w; }


        SIMD_INLINE IVector3D<T> GetV() const { return v; }
        SIMD_INLINE            T GetW() const { return w; }


        SIMD_INLINE T GetAngle() const { return T(2)*LogRotor().Length(); }


        /// Spinor operators : quantum mechanics
        SIMD_INLINE IQuaternion<T> GetConjugateSpinor()  const { return Quaternion(w,  x,  y,  z); }
        SIMD_INLINE IQuaternion<T> GetXConjugateSpinor() const { return Quaternion(w,  x, -y, -z); }
        SIMD_INLINE IQuaternion<T> GetYConjugateSpinor() const { return Quaternion(w, -x,  y, -z); }
        SIMD_INLINE IQuaternion<T> GetZConjugateSpinor() const { return Quaternion(w, -x, -y,  z); }



        /// Initialize the quaternion using Euler angles
        SIMD_INLINE void SetEulerAngles(T angleX, T angleY, T angleZ)
        {
            // Rotate angle Euler
            *this = FromEulerAngles(angleX,angleY,angleZ);

            // Normalize the quaternion
            Normalize();
        }


   SIMD_INLINE IVector3D<T> GetEulerAngles3() const
   {
        IVector3D<T> rpy;
        T a12, a22, a31, a32, a33;  // rotation matrix coefficients for Euler angles and gravity components
        a12 = 2.0f * (x * y + w * z);
        a22 = w * w + x * x - y * y - z * z;
        a31 = 2.0f * (w * x + y * z);
        a32 = 2.0f * (x * z - w * y);
        a33 = w * w - x * x - y * y + z * z;
        rpy[0] = atan2f(a31, a33);
        rpy[1] = -asinf(a32);
        rpy[2] = atan2f(a12, a22);
        rpy[0] *= 180.0f / PI;
        rpy[1] *= 180.0f / PI;
        rpy[2] *= 180.0f / PI;
        rpy[2] += 4.82f;

        if (rpy[2] >= +180.f)
            rpy[2] -= 360.f;
        else if (rpy[2] < -180.f)
            rpy[2] += 360.f;
            
        return rpy;
   }


 


        /// Quaternion using Euler angles
                SIMD_INLINE IVector3D<T> GetEulerAngles2() const
                {
                    // Euler angles in radians
                    T pitch;
                    T yaw;
                    T roll;

                    IQuaternion<T> q(*this);

                    T sqw = q.w * q.w;
                    T sqx = q.v.x * q.v.x;
                    T sqy = q.v.y * q.v.y;
                    T sqz = q.v.z * q.v.z;

                    // If quaternion is normalised the unit is one, otherwise it is the correction factor
                    T unit = sqx + sqy + sqz + sqw;
                    T val = q.v.x * q.v.y + q.v.z * q.w;

                    if (val > T(0.4999) * unit)                                // 0.4999f OR 0.5f - EPSILON
                    {
                        // Singularity at north pole
                        pitch = 2.f * IAtan2(q.v.x, q.w);          // Yaw
                        yaw = M_PI * 0.5f;                         // Pitch
                        roll = 0.f;                                // Roll
                    }
                    else if (val < -T(0.4999) * unit)                          // -0.4999f OR -0.5f + EPSILON
                    {
                        // Singularity at south pole
                        pitch = -2.f * IAtan2(q.v.x, q.w);         // Yaw
                        yaw = -M_PI * 0.5f;                        // Pitch
                        roll = 0.f;                                // Roll
                    }
                    else
                    {
                        pitch = IAtan2(2.f * q.v.y * q.w - 2.f * q.v.x * q.v.z,  sqx - sqy - sqz + sqw);   // Yaw
                        yaw   = IASin(2.f * val / unit);                                                   // Pitch
                        roll  = IAtan2(2.f * q.v.x * q.w - 2.f * q.v.y * q.v.z, -sqx + sqy - sqz + sqw);   // Roll
                    }

//                   pitch =  atan2(2*(v.y*v.z + w*v.x), w*w - v.x*v.x - v.y*v.y + v.z*v.z);
//                   yaw =  asin(-2*(v.x*v.z - w*v.y));
//                   roll = atan2(2*(v.x*v.y + w*v.z), w*w + v.x*v.x - v.y*v.y - v.z*v.z);

                    return /*NormalizeAngles*/(IVector3D<T>(roll,pitch,yaw));
                }




        /// Quaternion using Euler angles
        SIMD_INLINE IVector3D<T> GetEulerAngles() const
        {
            T pitch;
            T yaw;
            T roll;
            T val = 2.0f * (y * w - z * x);
            if (val  >= 0.99995f)
            {
                pitch = 2.0f * IAtan2(x, w);
                yaw = (0.5f  * M_PI);
                roll = 0.0f;
            }
            else if (val  <= -0.99995f)
            {
                pitch = 2.0f * IAtan2(x, w);
                yaw = -(0.5f * M_PI);
                roll = 0.0f;
            }
            else
            {
                yaw =   IASin (val);
                pitch = IAtan2 (2.0f * (x * w + z * y), 1.0f - 2.0f * (x * x + y * y));
                roll =  IAtan2 (2.0f * (z * w + x * y), 1.0f - 2.0f * (y * y + z * z));
            }
            return IVector3D<T>(pitch, yaw, roll);
        }



        ///////////////////////////////
        // Quaternion to Euler
        ///////////////////////////////
        enum RotSeq{zyx, zyz, zxy, zxz, yxz, yxy, yzx, yzy, xyz, xyx, xzy,xzx};


        // note:
        // return values of res[] depends on rotSeq.
        // i.e.
        // for rotSeq zyx,
        // x = res[0], y = res[1], z = res[2]
        // for rotSeq xyz
        // z = res[0], y = res[1], x = res[2]
        // ...
        IVector3D<T> GetEulerAngleGimbalLock(RotSeq rotSeq)
        {
            const IQuaternion<T> q(*this);
            T res[3];

            switch(rotSeq)
            {
            case zyx:
              threeaxisrot( 2*(q.x*q.y + q.w*q.z),
                             q.w*q.w + q.x*q.x - q.y*q.y - q.z*q.z,
                            -2*(q.x*q.z - q.w*q.y),
                             2*(q.y*q.z + q.w*q.x),
                             q.w*q.w - q.x*q.x - q.y*q.y + q.z*q.z,
                             res);
              break;

            case zyz:
              twoaxisrot( 2*(q.y*q.z - q.w*q.x),
                           2*(q.x*q.z + q.w*q.y),
                           q.w*q.w - q.x*q.x - q.y*q.y + q.z*q.z,
                           2*(q.y*q.z + q.w*q.x),
                          -2*(q.x*q.z - q.w*q.y),
                          res);
              break;

            case zxy:
              threeaxisrot( -2*(q.x*q.y - q.w*q.z),
                              q.w*q.w - q.x*q.x + q.y*q.y - q.z*q.z,
                              2*(q.y*q.z + q.w*q.x),
                             -2*(q.x*q.z - q.w*q.y),
                              q.w*q.w - q.x*q.x - q.y*q.y + q.z*q.z,
                              res);
              break;

            case zxz:
              twoaxisrot( 2*(q.x*q.z + q.w*q.y),
                          -2*(q.y*q.z - q.w*q.x),
                           q.w*q.w - q.x*q.x - q.y*q.y + q.z*q.z,
                           2*(q.x*q.z - q.w*q.y),
                           2*(q.y*q.z + q.w*q.x),
                           res);
              break;

            case yxz:
              threeaxisrot( 2*(q.x*q.z + q.w*q.y),
                             q.w*q.w - q.x*q.x - q.y*q.y + q.z*q.z,
                            -2*(q.y*q.z - q.w*q.x),
                             2*(q.x*q.y + q.w*q.z),
                             q.w*q.w - q.x*q.x + q.y*q.y - q.z*q.z,
                             res);
              break;

            case yxy:
              twoaxisrot( 2*(q.x*q.y - q.w*q.z),
                           2*(q.y*q.z + q.w*q.x),
                           q.w*q.w - q.x*q.x + q.y*q.y - q.z*q.z,
                           2*(q.x*q.y + q.w*q.z),
                          -2*(q.y*q.z - q.w*q.x),
                          res);
              break;

            case yzx:
              threeaxisrot( -2*(q.x*q.z - q.w*q.y),
                              q.w*q.w + q.x*q.x - q.y*q.y - q.z*q.z,
                              2*(q.x*q.y + q.w*q.z),
                             -2*(q.y*q.z - q.w*q.x),
                              q.w*q.w - q.x*q.x + q.y*q.y - q.z*q.z,
                              res);
              break;

            case yzy:
              twoaxisrot( 2*(q.y*q.z + q.w*q.x),
                          -2*(q.x*q.y - q.w*q.z),
                           q.w*q.w - q.x*q.x + q.y*q.y - q.z*q.z,
                           2*(q.y*q.z - q.w*q.x),
                           2*(q.x*q.y + q.w*q.z),
                           res);
              break;

            case xyz:
              threeaxisrot( -2*(q.y*q.z - q.w*q.x),
                            q.w*q.w - q.x*q.x - q.y*q.y + q.z*q.z,
                            2*(q.x*q.z + q.w*q.y),
                           -2*(q.x*q.y - q.w*q.z),
                            q.w*q.w + q.x*q.x - q.y*q.y - q.z*q.z,
                            res);
              break;

            case xyx:
              twoaxisrot( 2*(q.x*q.y + q.w*q.z),
                          -2*(q.x*q.z - q.w*q.y),
                           q.w*q.w + q.x*q.x - q.y*q.y - q.z*q.z,
                           2*(q.x*q.y - q.w*q.z),
                           2*(q.x*q.z + q.w*q.y),
                           res);
              break;

            case xzy:
              threeaxisrot( 2*(q.y*q.z + q.w*q.x),
                             q.w*q.w - q.x*q.x + q.y*q.y - q.z*q.z,
                            -2*(q.x*q.y - q.w*q.z),
                             2*(q.x*q.z + q.w*q.y),
                             q.w*q.w + q.x*q.x - q.y*q.y - q.z*q.z,
                             res);
              break;

            case xzx:
              twoaxisrot( 2*(q.x*q.z - q.w*q.y),
                           2*(q.x*q.y + q.w*q.z),
                           q.w*q.w + q.x*q.x - q.y*q.y - q.z*q.z,
                           2*(q.x*q.z + q.w*q.y),
                          -2*(q.x*q.y - q.w*q.z),
                          res);
              break;
            default:
              std::cout << "Unknown rotation sequence" << std::endl;
              break;
           }


            return IVector3D<T>(res[0],res[1],res[2]);
        }



        //        /// Quaternion using camera look direction
        //        SIMD_INLINE IQuaternion<T> LookRotation(IVector3D<T>& lookAt, IVector3D<T>& upDirection)
        //        {
        //            IVector3D<T> forward = lookAt.GetUnit();
        //            IVector3D<T> right   = upDirection.GetUnit().Cross(forward);
        //            IVector3D<T> up      = forward.Cross(right);

        //            IQuaternion<T> ret;
        //            ret.w = ISqrt(1.0f + right.x + up.y + forward.z) * 0.5f;
        //            T w4_recip = 1.0f / (4.0f * ret.w);
        //            ret.v.x = (forward.y - up.z) * w4_recip;
        //            ret.v.y = (right.z - forward.x) * w4_recip;
        //            ret.v.z = (up.x - right.y) * w4_recip;

        //            return ret;
        //        }

        /// Quaternion using camera look direction
        SIMD_INLINE IQuaternion<T> LookRotation(const IVector3D<T>& forward, const IVector3D<T>& up)
        {

            IVector3D<T> vector = (forward.Normalized());
            IVector3D<T> vector2 =(up.Cross(vector)).Normalized();
            IVector3D<T> vector3 =(vector.Cross(vector2));

            T m00 = vector2.x;
            T m01 = vector2.y;
            T m02 = vector2.z;
            T m10 = vector3.x;
            T m11 = vector3.y;
            T m12 = vector3.z;
            T m20 = vector.x;
            T m21 = vector.y;
            T m22 = vector.z;

            T num8 = (m00 + m11) + m22;
            IQuaternion<T> quaternion;
            if (num8 > T(0.f))
            {
                T num = ISqrt(num8 + T(1.f));
                quaternion.w = num * T(0.5f);
                num = T(0.5f) / num;
                quaternion.x = (m12 - m21) * num;
                quaternion.y = (m20 - m02) * num;
                quaternion.z = (m01 - m10) * num;
                return quaternion;
            }
            if ((m00 >= m11) && (m00 >= m22))
            {
                T num7 = ISqrt(((T(1.f) + m00) - m11) - m22);
                T num4 = T(0.5f) / num7;
                quaternion.x = T(0.5f) * num7;
                quaternion.y = (m01 + m10) * num4;
                quaternion.z = (m02 + m20) * num4;
                quaternion.w = (m12 - m21) * num4;
                return quaternion;
            }
            if (m11 > m22)
            {
                T num6 = ISqrt(((T(1.f) + m11) - m00) - m22);
                T num3 = T(0.5f) / num6;
                quaternion.x = (m10+ m01) * num3;
                quaternion.y = T(0.5f) * num6;
                quaternion.z = (m21 + m12) * num3;
                quaternion.w = (m20 - m02) * num3;
                return quaternion;
            }

            T num5 = ISqrt(((T(1.f) + m22) - m00) - m11);
            T num2 = T(0.5f) / num5;
            quaternion.x = (m20 + m02) * num2;
            quaternion.y = (m21 + m12) * num2;
            quaternion.z = T(0.5f) * num5;
            quaternion.w = (m01 - m10) * num2;
            return quaternion;
        }



        /// <summary>
        /// Creates a left-handed, look-at quaternion.
        /// </summary>
        /// <param name="eye">The position of the viewer's eye.</param>
        /// <param name="target">The camera look-at target.</param>
        /// <param name="up">The camera's up vector.</param>
        /// <param name="result">When the method completes, contains the created look-at quaternion.</param>
        static const IQuaternion<T> LookAtLH(const IVector3D<T>& eye, const IVector3D<T>& target, const IVector3D<T>& up)
        {
            return IQuaternion<T>(IMatrix3x3<T>::LookAtLH(eye, target, up));
        }


        /// <summary>
        /// Creates a right-handed, look-at quaternion.
        /// </summary>
        /// <param name="eye">The position of the viewer's eye.</param>
        /// <param name="target">The camera look-at target.</param>
        /// <param name="up">The camera's up vector.</param>
        /// <param name="result">When the method completes, contains the created look-at quaternion.</param>
        static const IQuaternion<T> LookAtRH(const IVector3D<T>& eye, const IVector3D<T>& target, const IVector3D<T>& up)
        {
            return IQuaternion<T>(IMatrix3x3<T>::LookAtRH(eye, target, up));
        }



        //--------------------- operators -------------------------//

        /**
         * Copy operator
         * @param rhs Right hand side argument of binary operator.
         */
        SIMD_INLINE IQuaternion<T>& operator=(const IQuaternion<T>& rhs)
        {
            x = rhs.x;
            y = rhs.y;
            z = rhs.z;
            w = rhs.w;
            return *this;
        }

        /**
         * Copy convert operator
         * @param rhs Right hand side argument of binary operator.
         */
        template<class FromT>
        SIMD_INLINE IQuaternion<T>& operator=(const IQuaternion<FromT>& rhs)
        {
            x = static_cast<T>(rhs.x);
            y = static_cast<T>(rhs.y);
            z = static_cast<T>(rhs.z);
            w = static_cast<T>(rhs.w);
            return *this;
        }


        /**
         * Addition operator
         * @param rhs Right hand side argument of binary operator.
         */
        SIMD_INLINE IQuaternion<T> operator+(const IQuaternion<T>& rhs) const
        {
            const IQuaternion<T>& lhs = *this;
            return IQuaternion<T>(lhs.w + rhs.w,
                                  lhs.x + rhs.x,
                                  lhs.y + rhs.y,
                                  lhs.z + rhs.z );
        }


        /**
         * Subtraction operator
         * @param rhs Right hand side argument of binary operator.
         */
        SIMD_INLINE IQuaternion<T> operator-(const IQuaternion<T>& rhs) const
        {
            const IQuaternion<T>& lhs = *this;
            return IQuaternion<T>(lhs.w - rhs.w,
                                  lhs.x - rhs.x,
                                  lhs.y - rhs.y,
                                  lhs.z - rhs.z);
        }


        /**
         * Multiplication operator
         * @param rhs Right hand side argument of binary operator.
         */
        SIMD_INLINE IQuaternion<T> operator*(T rhs) const
        {
            return IQuaternion<T>(w * rhs,
                                  x * rhs,
                                  y * rhs,
                                  z * rhs);
        }


        /**
         * Multiplication invert operator
         * @param rhs Right hand side argument of binary operator.
         */
        SIMD_INLINE IQuaternion<T> operator/(T rhs) const
        {
            return IQuaternion<T>(w / rhs,
                                  x / rhs,
                                  y / rhs,
                                  z / rhs);
        }



        /**
          * Quaternion multiplication operator
          * @param rhs Right hand side argument of binary operator.
          */
        //        SIMD_INLINE IQuaternion<T> operator*(const IQuaternion<T>& rhs) const
        //        {
        //            return IQuaternion<T> (w * rhs.w - v.Dot(rhs.v), w * rhs.v + rhs.w * v + v.Cross(rhs.v));
        //        }


        /**
        * Quaternion multiplication operator
        * @param rhs Right hand side argument of binary operator.
        */
        SIMD_INLINE IQuaternion<T> operator*(const IQuaternion<T>& Q) const
        {
            return IQuaternion(w*Q.w - x*Q.x - y*Q.y - z*Q.z,
                               w*Q.x + x*Q.w + y*Q.z - z*Q.y,
                               w*Q.y - x*Q.z + y*Q.w + z*Q.x,
                               w*Q.z + x*Q.y - y*Q.x + z*Q.w);
        }



        /**
        * Quaternion multiplication operator
        * @param rhs Right hand side argument of binary operator.
        */
        SIMD_INLINE IVector3D<T> operator*(const IVector3D<T>& rhs) const
        {
            IVector3D<T> u(x, y, z);
            return u * (u.Dot(rhs) * T(2.0)) + rhs * (w * w - u.Dot(u)) + u.Cross(rhs) * (T(2.0) * w);
        }



        /// multiplication operator dot vector
        SIMD_INLINE IVector3D<T> nVidia_dot(const IVector3D<T>& _v) const
        {
            // nVidia SDK implementation
            IVector3D<T> uv, uuv;
            IVector3D<T> qvec(x, y, z);
            uv = qvec.Cross(_v);
            uuv = qvec.Cross(uv);
            uv *= (2.0f * w);
            uuv *= 2.0f;
            return _v + uv + uuv;
        }


        //----------------------------------------------------------//

        /**
        * Addition operator
        * @param rhs Right hand side argument of binary operator.
        */
        SIMD_INLINE IQuaternion<T>& operator+=(const IQuaternion<T>& rhs)
        {
            w += rhs.w;
            x += rhs.x;
            y += rhs.y;
            z += rhs.z;
            return *this;
        }

        /**
        * Subtraction operator
        * @param rhs Right hand side argument of binary operator.
        */
        SIMD_INLINE IQuaternion<T>& operator-=(const IQuaternion<T>& rhs)
        {
            w -= rhs.w;
            x -= rhs.x;
            y -= rhs.y;
            z -= rhs.z;
            return *this;
        }


        /**
        * Multiplication operator
        * @param rhs Right hand side argument of binary operator.
        */
        SIMD_INLINE IQuaternion<T>& operator*=(T rhs)
        {
            w *= rhs;
            x *= rhs;
            y *= rhs;
            z *= rhs;
            return *this;
        }


        /**
        * Multiplication invert operator
        * @param rhs Right hand side argument of binary operator.
        */
        SIMD_INLINE IQuaternion<T>& operator/=(T rhs)
        {
            w /= rhs;
            x /= rhs;
            y /= rhs;
            z /= rhs;
            return *this;
        }


        /**
        * Multiplication operator
        * @param rhs Right hand side argument of binary operator.
        */
        SIMD_INLINE IQuaternion<T>& operator*=(const IQuaternion<T>& rhs)
        {
            IQuaternion q = (*this) * rhs;
            w = q.w;
            x = q.x;
            y = q.y;
            z = q.z;
            return *this;
        }



        /**
        * Equality test operator
        * @param rhs Right hand side argument of binary operator.
        * @note Test of equality is based of threshold EPSILON value. To be two
        * values equal, must satisfy this condition | lhs - rhs | < EPSILON,
        * for all IQuaternion coordinates.
        */
        SIMD_INLINE bool operator==(const IQuaternion<T>& rhs) const
        {
            const IQuaternion<T>& lhs = *this;
            return (IAbs(lhs.w - rhs.w) < MACHINE_EPSILON) &&
                   (IAbs(lhs.x - rhs.x) < MACHINE_EPSILON) &&
                   (IAbs(lhs.y - rhs.y) < MACHINE_EPSILON) &&
                   (IAbs(lhs.z - rhs.z) < MACHINE_EPSILON);
        }


        /**
        * Inequality test operator
        * @param rhs Right hand side argument of binary operator.
        * @return not (lhs == rhs) :-P
        */
        SIMD_INLINE bool operator!=(const IQuaternion<T>& rhs) const
        {
            return !(*this == rhs);
        }


        //---------------------- Methods ---------------------//


        SIMD_INLINE IQuaternion<T> Cross(const IQuaternion<T>& Q) const
        {
            return IQuaternion<T>(0, -z*Q.y+y*Q.z,
                                      z*Q.x-x*Q.z,
                                     -y*Q.x+x*Q.y);
        }

        SIMD_INLINE IQuaternion<T> Dot(const IQuaternion<T>& Q) const
        {
            return x*Q.x + y*Q.y + z*Q.z;
        }


        SIMD_INLINE IQuaternion<T> Commutator(const IQuaternion<T>& Q) const
        {
            return IQuaternion<T>(0, -2*z*Q.y+2*y*Q.z,
                                      2*z*Q.x-2*x*Q.z,
                                     -2*y*Q.x+2*x*Q.y);
        }


        SIMD_INLINE IQuaternion<T> Pow(const T t) const
        {
            return (this->Log() * t).Exp();
        }

        SIMD_INLINE IQuaternion<T> Pow(const IQuaternion<T>& Q) const
        {
            return (this->Log() * Q).Exp();
        }

     //-------------[ unary operations ]--------------------------

     /**
     * Unary negate operator
     * @return negated IQuaternion
     */
     SIMD_INLINE IQuaternion<T> operator-() const
     {
        return IQuaternion<T>( -w, -x , -y, -z );
     }

     /**
     * Unary conjugate operator
     * @return conjugated IQuaternion
     */
     SIMD_INLINE IQuaternion<T> operator~() const
     {
        return IQuaternion<T>( w, -x, -y, -z );
     }

     /**
     * Get lenght of IQuaternion.
     * @return Length of IQuaternion.
     */
     SIMD_INLINE T Length() const
     {
        return ISqrt(w * w + x * x + y * y + z * z);
     }

     /**
     * Return square of length.
     * @return length ^ 2
     * @note This method is faster then length(). For comparison
     * of length of two IQuaternion can be used just this value, instead
     * of more expensive length() method.
     */
     SIMD_INLINE T LengthSquare() const
     {
        return w * w + x * x + y * y + z * z;
     }

     /**
     * Normalize IQuaternion
     */
     SIMD_INLINE void Normalize()
     {
         T len = Length();
         w /= len;
         x /= len;
         y /= len;
         z /= len;
     }

     /**
     * Inverse IQuaternion
     */
     SIMD_INLINE IQuaternion<T> GetInverse() const
     {
         T lengthSquareQuaternion = LengthSquare();

         assert (lengthSquareQuaternion > MACHINE_EPSILON);

         // Compute and return the inverse quaternion
         return IQuaternion<T>( w / lengthSquareQuaternion,
                               -x / lengthSquareQuaternion,
                               -y / lengthSquareQuaternion,
                               -z / lengthSquareQuaternion);

     }


     /**
     * Unit IQuaternion
     */
     SIMD_INLINE IQuaternion<T> GetUnit() const
     {
         T lengthQuaternion = Length();

         // Check if the length is not equal to zero
         assert (lengthQuaternion > MACHINE_EPSILON);

         // Compute and return the unit quaternion
         return IQuaternion<T>(w / lengthQuaternion,
                               x / lengthQuaternion,
                               y / lengthQuaternion,
                               z / lengthQuaternion);
     }


    /**
    * Conjugate IQuaternion
    */
    SIMD_INLINE IQuaternion<T> GetConjugate() const
    {
       return IQuaternion<T>(w , -x , -y , -z);
    }



    /**
        * Creates IQuaternion for eulers angles.
        * @param x Rotation around x axis (in degrees).
        * @param y Rotation around y axis (in degrees).
        * @param z Rotation around z axis (in degrees).
        * @return IQuaternion object representing transformation.
        */
        static SIMD_INLINE IQuaternion<T> FromEulerAngles(T x, T y, T z)
        {
               IQuaternion<T> ret = FromAngleAxis(IVector3D<T>(1, 0, 0), x) *
                                    FromAngleAxis(IVector3D<T>(0, 1, 0), y) *
                                    FromAngleAxis(IVector3D<T>(0, 0, 1), z);
               return ret;
        }



        /**
        * Creates IQuaternion as rotation around axis.
        * @param axis Unit vector expressing axis of rotation.
        * @param angleDeg Angle of rotation around axis (in degrees).
        */
        static SIMD_INLINE IQuaternion<T> FromAngleAxis(const IVector3D<T> axis, T radianDegrees)
        {
               T angleRad = /*IDegreesToRadians*/(radianDegrees) * T(0.5);
               T sa2 = ISin(angleRad);
               T ca2 = ICos(angleRad);
               return IQuaternion<T>( ca2 , sa2 * axis );
        }



     /**
     * Creates IQuaternion as rotation around axis.
     * @param axis Unit vector expressing axis of rotation.
     * @param angleDeg Angle of rotation around axis (in degrees).
     */
     static SIMD_INLINE IQuaternion<T> FromAxisRot(const IVector3D<T> axis, T angleDeg)
     {
            T angleRad = IDegreesToRadians(angleDeg) * T(0.5);
            T sa2 = ISin(angleRad);
            T ca2 = ICos(angleRad);
            return IQuaternion<T>( ca2 , sa2 * axis );
     }


    //        IQuaternion<T>  createRotation( IVector3D<T> RotAxis , T AngleInRads )
    //        {
    //            T Thetad2 = AngleInRads * T(0.5);
    //            T SinTd2 = ISin(Thetad2);
    //            T CosTd2 = ICos(Thetad2);

    //            // Normalize the quaternion axis
    //            RotAxis.Normalized();

    //            w = CosTd2;
    //            x = (SinTd2 * RotAxis).x;
    //            y = (SinTd2 * RotAxis).y;
    //            z = (SinTd2 * RotAxis).z;

    //            return *this;
    //        }








        /**
        * Linear interpolation of two IQuaternions
        * @param fact Factor of interpolation. For translation from position
        * of this vector to IQuaternion rhs, values of factor goes from 0.0 to 1.0.
        * @param rhs Second IQuaternion for interpolation
        * @note However values of fact parameter are reasonable only in interval
        * [0.0 , 1.0], you can pass also values outside of this interval and you
        * can Get result (extrapolation?)
        */
        SIMD_INLINE IQuaternion<T> Lerp(T fact, const IQuaternion<T>& rhs) const
        {
             return IQuaternion<T>((1 - fact) * w + fact * rhs.w, v.Lerp(fact, rhs.v));
        }


        /**
        * Computes spherical interpolation between IQuaternions (this, q2)
        * using coefficient of interpolation r (in [0, 1]).
        *
        * @param r The ratio of interpolation form this (r = 0) to q2 (r = 1).
        * @param q2 Second IQuaternion for interpolation.
        * @return Result of interpolation.
        */
        SIMD_INLINE IQuaternion<T> Slerp(T r, const IQuaternion<T>& q2) const
        {
            IQuaternion<T> ret;
            T cosTheta = w * q2.w + x * q2.x + y * q2.y + z * q2.z;
            T theta = IACos(cosTheta);
            if (IAbs(theta) < MACHINE_EPSILON)
            {
                ret = *this;
            }
            else
            {
                T sinTheta = ISqrt(1.0 - cosTheta * cosTheta);
                if (IAbs(sinTheta) < MACHINE_EPSILON)
                {
                    ret.w = 0.5 * w + 0.5 * q2.w;
                    ret.v = GetImage().Lerp(0.5, q2.GetImage());
                }
                else
                {
                    T rA = ISin((1.0 - r) * theta) / sinTheta;
                    T rB = ISin(r * theta) / sinTheta;

                    ret.w = w * rA + q2.w * rB;
                    ret.x = x * rA + q2.x * rB;
                    ret.y = y * rA + q2.y * rB;
                    ret.z = z * rA + q2.z * rB;
                }
            }
            return ret;
        }


        /// The angle between the vectors is simple to find: the dot product gives its cosine.
        /// The needed axis is also simple to find: it’s the cross product of the two vectors.
        static SIMD_INLINE IQuaternion<T> SlerpToVectors(  IVector3D<T> start ,  IVector3D<T> dest )
        {

            //start.normalize();
            //dest.normalize();

            T cosTheta = start.Dot(dest);
            IVector3D<T> rotationAxis;

            if (cosTheta < -1 + 0.0001f)
            {
                // special case when vectors in opposite directions:
                // there is no "ideal" rotation axis
                // So guess one; any will do as long as it's perpendicular to start
                rotationAxis = IVector3D<T>::Z.Cross(start);
                if (rotationAxis.LengthSquare() < T(0.0001) )
                {
                    // bad luck, they were parallel, try again!
                    rotationAxis = IVector3D<T>::X.Cross(start);
                }

                rotationAxis.Normalize();
                static IQuaternion<T> q;
                return q.FromAxisRot(rotationAxis,T(180.0f));
            }

            rotationAxis = start.Cross(dest);
            rotationAxis.Normalize();


            T s = ISqrt( (1+cosTheta)*2 );
            T invs = 1.0 / s;

            return IQuaternion<T> ( s * 0.5f , rotationAxis * invs);

        }




        /// Return logarithm of Quaternion.
        SIMD_INLINE IQuaternion<T> Log() const
        {
            IQuaternion<T>  Result;
            const T b = ISqrt(x*x + y*y + z*z);
            if(IAbs(b) <= MACHINE_EPSILON*IAbs(w))
            {
                if(w<0.0)
                {
                    //INFOTOCERR << " Error: Infinitely many solutions for log of a negative scalar (w=" << w << ")." << endl;
                    //throw(InfinitelyManySolutions);
                }
                Result.w = ILog(w);
            }
            else
            {
                const T _v = IAtan2(b, w);
                const T _f = _v/b;
                Result.w   = ILog(w*w+b*b)/2.0;
                // Result.w = std::log(w/std::cos(v)); // Not nice for unit vectors [w=cos(v)=0]
                Result.x = _f*x;
                Result.y = _f*y;
                Result.z = _f*z;
            }
            return Result;
        }



        /// Return logarithm of a rotor.
        SIMD_INLINE IQuaternion<T> LogRotor() const
        {
            /// This function is just like the standard log function, except
            /// that a negative scalar does not raise an exception; instead, the
            /// scalar pi is returned.
            IQuaternion<T> Result;
            const T b = ISqrt(x*x + y*y + z*z);

            if(IAbs(b) <= MACHINE_EPSILON*IAbs(w))
            {
                if(w<0.0)
                {
                    //    INFOTOCERR
                    //   << "\nWarning: Infinitely many solutions for log of a negative scalar (w=" << w << "); "
                    //   << "arbitrarily returning the one in the x direction." << endl;
                    Result.x = M_PI;
                    if(IAbs(w+1)>MACHINE_EPSILON)
                    {
                        Result.w = ILog(-w);
                    }
                }
                else
                {
                    Result.w = ILog(w);
                }
            }
            else
            {
                const T _v = IAtan2(b, w);
                const T _f = _v/b;
                Result.w = ILog(w*w+b*b)/2.0;
                // Result.w = std::log(w/std::cos(v)); // Not nice for unit vectors [w=cos(v)=0]
                Result.x = _f*x;
                Result.y = _f*y;
                Result.z = _f*z;
            }

            return Result;
        }

        /// Return exponent of Quaternion.
        SIMD_INLINE IQuaternion<T>  Exp() const
        {
            IQuaternion<T>  Result;
            const T b = ISqrt(x*x + y*y + z*z);
            if(IAbs(b)<=MACHINE_EPSILON*IAbs(w))
            {
                Result.w = IExp(w);
            }
            else
            {
                const T e = IExp(w);
                const T f = ISin(b)/b; // Note: b is never 0.0 at this point
                Result.w   = e*ICos(b);
                Result.x = e*f*x;
                Result.y = e*f*y;
                Result.z = e*f*z;
            }
            return Result;
        }



        //----------[ output operator ]----------------------------
        /**
        * Provides output to standard output stream.
        */
        friend std::ostream& operator <<(std::ostream& oss, const IQuaternion<T>& q)
        {
            oss << "(" << "Re: " << q.w << " Im: " << q.v << ")";
            return oss;
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


        //-----------------------------------------------------------------------//


        //------------------------------- friendship ----------------------------//
        friend class IMatrix3x3<T>;
        friend class IMatrix4x4<T>;



    public:


    /**
    * The multiplicitive identity quaternion
    */
    static const IQuaternion<T> IDENTITY;

    /**
     * The additive identity quaternion.
     */
     static const IQuaternion<T> ZERO;

};


template<class T> const IQuaternion<T> IQuaternion<T>::IDENTITY(1.0, 0.0, 0.0, 0.0);
template<class T> const IQuaternion<T> IQuaternion<T>::ZERO(0.0, 0.0, 0.0, 0.0);



//--------------------------------------
// Typedef shortcuts for Quaternion
//-------------------------------------


using IQuaternionr    = IQuaternion<Real>;
using IQuaternionf    = IQuaternion<float>;
using IQuaterniond    = IQuaternion<double>;
using IQuaternioni    = IQuaternion<std::int32_t>;
using IQuaternionui   = IQuaternion<std::uint32_t>;
using IQuaternionb    = IQuaternion<std::int8_t>;
using IQuaternionub   = IQuaternion<std::uint8_t>;

} /* namespace  */



#endif /* IQUATERNION_H_ */
