 /********************************************************************************
 *
 * IPlane.h
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


#ifndef IPLANE_H
#define IPLANE_H

#include "IVector3D.h"
#include "IQuaternion.h"


namespace IMath
{

//-------------------------------------------------------------------------------
//-- Classes --------------------------------------------------------------------
//-------------------------------------------------------------------------------
template<class T>
class IPlane
{

protected:


   // IVector3D<T>  mOrigin;
    IVector3D<T>  mNormal;
    T             mOffset;


public:
    // constructor/destructor

    // Set this plane from a normal and a signed distance from its origin.
    IPlane( const IVector3D<T>& normal = IVector3D<T>(1,0,0) , T offset = T(0) )
      : mNormal(normal) ,
        mOffset(offset)
    {

    }

    // Set this plane from a normal and a point on the plane.
    IPlane(const IVector3D<T>& _normal, const IVector3D<T>& _point)
    {
        mNormal = _normal;
        mOffset = _normal.Dot(_point);
       // mOrigin = _point;
    }

    IPlane( const IVector3D<T>& p0,
            const IVector3D<T>& p1,
            const IVector3D<T>& p2 )
    {
       Set( p0, p1, p2 );
    }


    IPlane( T a, T b, T c, T d )
    {
       Set( a, b, c, d );
    }

    // copy operations
    IPlane(const IPlane& other)
     : mNormal( other.mNormal ),
       mOffset( other.mOffset )
    {

    }


    SIMD_INLINE ~IPlane() {}


    // ---------------------------------------------------------------------------
    // Assigment operator
    //-----------------------------------------------------------------------------
    IPlane& operator=(const IPlane& other)
    {
        // if same object
        if ( this == &other )
            return *this;

        mNormal = other.mNormal;
        mOffset = other.mOffset;

        return *this;
    }


    // accessors
    SIMD_INLINE const IVector3D<T>& GetNormal() const { return mNormal; }
    SIMD_INLINE T GetOffset() const { return mOffset; }

    // ---------------------------------------------------------------------------
    // Returns the two endpoints
    //-----------------------------------------------------------------------------
    void Get( IVector3D<T>& normal, T& offset ) const
    {
        normal = mNormal;
        offset = mOffset;
    }

    // ---------------------------------------------------------------------------
    // Are two IPlane's equal?
    //----------------------------------------------------------------------------
    bool operator==( const IPlane&  plane ) const
    {
        return (plane.mNormal == mNormal &&
                plane.mOffset == mOffset);
    }

    // ---------------------------------------------------------------------------
    // Are two IPlane's not equal?
    //----------------------------------------------------------------------------
    bool operator!=( const IPlane&  plane ) const
    {
        return !(plane.mNormal == mNormal &&
                 plane.mOffset == mOffset);
    }

    // manipulators
    SIMD_INLINE void Set( const IVector3D<T>& n, T d )
    {
        Set( n.x, n.y, n.z, d );
    }

    // ---------------------------------------------------------------------------
    // Sets the parameters
    //-----------------------------------------------------------------------------
    void Set( T a, T b, T c, T d )
    {
        // normalize for cheap distance checks
        T lensq = a*a + b*b + c*c;
        // length of normal had better not be zero
        assert( !iIsZero( lensq ) );

        // recover gracefully
        if ( iIsZero( lensq ) )
        {
            mNormal = IVector3D<T>::X;
            mOffset = 0.0f;
        }
        else
        {
            T recip = 1.0/iSqrt(lensq);
            mNormal.setAllValues( a*recip, b*recip, c*recip );
            mOffset = d*recip;
        }
    }



    // ---------------------------------------------------------------------------
    // Sets the parameters
    //-----------------------------------------------------------------------------
    void Set( const IVector3D<T>& p0, const IVector3D<T>& p1, const IVector3D<T>& p2 )
    {
        mNormal = (IVector3D<T>::triNormal(p0,p1,p2));
        mOffset = mNormal.Dot(p0);
    }

    // ---------------------------------------------------------------------------
    // Transforms plane into new space
    //-----------------------------------------------------------------------------
    IPlane Transform( T scale, const IQuaternion<T>& rotate, const IVector3D<T>& translate ) const
    {
        IPlane<T> plane;

        // get rotation matrix
        IMatrix3x3<T>  rotmatrix = rotate.getMatrix();

        // transform to get normal
        plane.mNormal = rotmatrix*mNormal/scale;

        // transform to get offset
        IVector3D<T> newTrans = translate*rotmatrix;
        plane.mOffset = -newTrans.Dot( mNormal )/scale + mOffset;

        return plane;
    }

    // ---------------------------------------------------------------------------
    // Transforms plane into new space
    //-----------------------------------------------------------------------------
    IPlane Transform( T scale, const IMatrix3x3<T>& rotmatrix, const IVector3D<T>& translate ) const
    {
        IPlane<T> plane;

        // transform to get normal
        plane.mNormal = rotmatrix*mNormal/scale;

        // transform to get offset
        IVector3D<T> newTrans = translate*rotmatrix;
        plane.mOffset = -newTrans.Dot( mNormal )/scale + mOffset;

        return plane;
    }

    // distance
    static T Distance( const IPlane& plane, const IVector3D<T>& point )
    {
        return ( plane.Test( point ) );
    }

    // ---------------------------------------------------------------------------
    // Returns the closest point on plane to point
    //-----------------------------------------------------------------------------
    IVector3D<T> ClosestPoint( const IVector3D<T>& point ) const
    {
        return point - Test(point)*mNormal;
    }

    // result of plane test
    SIMD_INLINE T Test( const IVector3D<T>& point ) const
    {
        return mNormal.Dot(point) - mOffset;
    }

    // result of plane test
    SIMD_INLINE T InvTest( const IVector3D<T>& point ) const
    {
        return mOffset - mNormal.Dot(point);
    }




    // vIntersectionLineToRay(): find the 3D intersection of a segment and a plane
    //    Input:  S = a segment, and Pn = a plane = {Point V0;  Vector n;}
    //    Output: *I0 = the intersect point (when it exists)
    IVector3D<T> VIntersectionRayToPlane( const IVector3D<T>& _ray_origin , const IVector3D<T>& _ray_dir ) const
    {
        IVector3D<T> N = mNormal;
        IVector3D<T> P = _ray_origin;
        IVector3D<T> W = _ray_dir;

        T  d =  InvTest(P);
        T  e =  N.Dot(W);

        if( IAbs(e) < MACHINE_EPSILON  ) return P;

        T param = d/e;

        return P + W * param;
    }



    //----------[ output operator ]----------------------------
     /**
     * Provides output to standard output stream.
     */
     friend std::ostream& operator <<(std::ostream& oss, const IPlane<T>& rhs)
     {
         oss << "(" << "normal: " << rhs.mNormal << " offset: " << rhs.mOffset << ")";
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


};



//--------------------------------------
// Typedef shortcuts for IPlane
//-------------------------------------

using IPlaner    = IPlane<Real>;
using IPlanef    = IPlane<float>;
using IPlaned    = IPlane<double>;
using IPlanei    = IPlane<std::int32_t>;
using IPlaneui   = IPlane<std::uint32_t>;
using IPlaneb    = IPlane<std::int8_t>;
using IPlaneub   = IPlane<std::uint8_t>;

}


#endif // IPLANE_H
