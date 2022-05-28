
 /********************************************************************************
 *
 * IOctonion.h
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

#ifndef IOCTONION_HPP_
#define IOCTONION_HPP_


#include "IComplex.h"
#include "IQuaternion.h"

namespace IMath
{

/**
 * Algebra Kelly
 */
template<class T> class IOctonion
{
	private:

		T   a;
		T   b;
		T   c;
		T   d;
		T   e;
		T   f;
		T   g;
		T   h;

	public:

   SIMD_INLINE IOctonion(T const& requested_a ,
						  T const& requested_b ,
						  T const& requested_c ,
						  T const& requested_d ,
						  T const& requested_e ,
						  T const& requested_f ,
						  T const& requested_g ,
						  T const& requested_h )
	 : a(requested_a),
	   b(requested_b),
	   c(requested_c),
	   d(requested_d),
	   e(requested_e),
	   f(requested_f),
	   g(requested_g),
	   h(requested_h)
	{
				// nothing to do!
	}

	 // constructor for H seen as C^4
   SIMD_INLINE IOctonion( IComplex<T> z0,
						   IComplex<T> z1 = IComplex<T>(),
						   IComplex<T> z2 = IComplex<T>(),
						   IComplex<T> z3 = IComplex<T>())
     :   a(z0.GetReal()),
         b(z0.GetImag()),
         c(z1.GetReal()),
         d(z1.GetImag()),
         e(z2.GetReal()),
         f(z2.GetImag()),
         g(z3.GetReal()),
         h(z3.GetImag())
	 {
		 // nothing to do!
	 }


	 // constructor for O seen as H^2
   SIMD_INLINE IOctonion( IQuaternion<T> q0,
                          IQuaternion<T> q1 = IQuaternion<T>())
	 :   a(q0.getW()),
	     b(q0.getV().x),
	     c(q0.getV().y),
	     d(q0.getV().z),
	     e(q1.getW()),
	     f(q1.getV().x),
	     g(q1.getV().y),
	     h(q1.getV().z)
	 {
		 // nothing to do!
	 }


	 template<typename X>
	 SIMD_INLINE  IOctonion(IOctonion<X> a_recopier)
	 :   a(static_cast<T>(a_recopier.R_component_1())),
		 b(static_cast<T>(a_recopier.R_component_2())),
		 c(static_cast<T>(a_recopier.R_component_3())),
		 d(static_cast<T>(a_recopier.R_component_4())),
		 e(static_cast<T>(a_recopier.R_component_5())),
		 f(static_cast<T>(a_recopier.R_component_6())),
		 g(static_cast<T>(a_recopier.R_component_7())),
		 h(static_cast<T>(a_recopier.R_component_8()))
	 {
		 // nothing to do!
	 }


     SIMD_INLINE T            GetReal()   const { return(a); }
     SIMD_INLINE IOctonion<T> GetUnreal() const { return( IOctonion<T>(static_cast<T>(0),b,c,d,e,f,g,h)); }

	 // number theory //
	 SIMD_INLINE T R_component_1() const { return(a); }
	 SIMD_INLINE T R_component_2() const { return(b); }
	 SIMD_INLINE T R_component_3() const { return(c); }
	 SIMD_INLINE T R_component_4() const { return(d); }
	 SIMD_INLINE T R_component_5() const { return(e); }
	 SIMD_INLINE T R_component_6() const { return(f); }
	 SIMD_INLINE T R_component_7() const { return(g); }
	 SIMD_INLINE T R_component_8() const { return(h); }

	 SIMD_INLINE IComplex<T> C_component_1() const { return(IComplex<T>(a,b)); }
	 SIMD_INLINE IComplex<T> C_component_2() const { return(IComplex<T>(c,d)); }
	 SIMD_INLINE IComplex<T> C_component_3() const { return(IComplex<T>(e,f)); }
	 SIMD_INLINE IComplex<T> C_component_4() const { return(IComplex<T>(g,h)); }

	 SIMD_INLINE IQuaternion<T> H_component_1() const { return(IQuaternion<T>(a,b,c,d)); }
	 SIMD_INLINE IQuaternion<T> H_component_2() const { return(IQuaternion<T>(e,f,g,h)); }


	 //---------------- operators --------------------//

	 // other assignment-related operators
	 //
	 // NOTE:    Octonion multiplication is *NOT* commutative;
	 //            symbolically, "q *= rhs;" means "q = q * rhs;"
	 //            and "q /= rhs;" means "q = q * inverse_of(rhs);";
	 //            octonion multiplication is also *NOT* associative


	 //--------------------------------------------------//

	 SIMD_INLINE IOctonion<T>&  operator += (const T& rhs)
	 {
         T    at = a + rhs; // exception guard
              a = at;
		 return(*this);
	 }

	 SIMD_INLINE IOctonion<T>  operator + (const T& rhs) const
	 {
		 T    at = a + rhs;    // exception guard
		 return IOctonion<T>( at , b , c , d ,e ,f ,g ,h);
	 }


	 //--------------------------------------------------//

	 SIMD_INLINE IOctonion<T>&   operator += (const IComplex<T>& rhs)
	 {
         T    at = a + rhs.GetReal();    // exception guard
         T    bt = b + rhs.GetImag();    // exception guard

		 a = at;
		 b = bt;

		 return(*this);
	 }

	 SIMD_INLINE IOctonion<T>   operator + (const IComplex<T>& rhs) const
	 {
         T    at = a + rhs.GetReal();    // exception guard
         T    bt = b + rhs.GetImag();    // exception guard

		 return IOctonion<T>(at,bt,c,d,e,f,g,h);
	 }

	 //--------------------------------------------------//

	 SIMD_INLINE IOctonion<T>&  operator += (const IQuaternion<T>& rhs)
	 {
		 T    at = a + rhs.R_component_1();    // exception guard
		 T    bt = b + rhs.R_component_2();    // exception guard
		 T    ct = c + rhs.R_component_3();    // exception guard
		 T    dt = d + rhs.R_component_4();    // exception guard

		 a = at;
		 b = bt;
		 c = ct;
		 d = dt;

		 return(*this);
	 }

	 SIMD_INLINE IOctonion<T>  operator + (const IQuaternion<T>& rhs)
	 {
		 T    at = a + rhs.R_component_1();    // exception guard
		 T    bt = b + rhs.R_component_2();    // exception guard
		 T    ct = c + rhs.R_component_3();    // exception guard
		 T    dt = d + rhs.R_component_4();    // exception guard

		 return IOctonion<T>(at,bt,ct,dt,e,f,g,h);
	 }

	 //--------------------------------------------------//


	 SIMD_INLINE IOctonion<T>&  operator += (const IOctonion<T>& rhs)
	 {
		 T    at = a + (rhs.R_component_1());    // exception guard
		 T    bt = b + (rhs.R_component_2());    // exception guard
		 T    ct = c + (rhs.R_component_3());    // exception guard
		 T    dt = d + (rhs.R_component_4());    // exception guard
		 T    et = e + (rhs.R_component_5());    // exception guard
		 T    ft = f + (rhs.R_component_6());    // exception guard
		 T    gt = g + (rhs.R_component_7());    // exception guard
		 T    ht = h + (rhs.R_component_8());    // exception guard

		 a = at;
		 b = bt;
		 c = ct;
		 d = dt;
		 e = et;
		 f = ft;
		 g = gt;
		 h = ht;

		 return(*this);
	 }

	 SIMD_INLINE IOctonion<T>  operator + (const IOctonion<T>& rhs) const
	 {
		 T    at = a + (rhs.R_component_1());    // exception guard
		 T    bt = b + (rhs.R_component_2());    // exception guard
		 T    ct = c + (rhs.R_component_3());    // exception guard
		 T    dt = d + (rhs.R_component_4());    // exception guard
		 T    et = e + (rhs.R_component_5());    // exception guard
		 T    ft = f + (rhs.R_component_6());    // exception guard
		 T    gt = g + (rhs.R_component_7());    // exception guard
		 T    ht = h + (rhs.R_component_8());    // exception guard

		 IOctonion<T> oc(at,bt,ct,dt,et,ft,gt,ht);
		 return oc;
	 }

	 //--------------------------------------------------//

	 SIMD_INLINE IOctonion<T>& operator -= (const T& rhs)
	 {
		 T    at = a - rhs;    // exception guard
		 a = at;

		 return(*this);
	 }

	 SIMD_INLINE IOctonion<T>  operator - (const T& rhs) const
	 {
		 T    at = a - rhs;    // exception guard

		 return IOctonion<T>( at , b , c , d ,e ,f ,g ,h);
	 }

	 //--------------------------------------------------//

	 SIMD_INLINE IOctonion<T>&  operator -= (const IComplex<T>& rhs)
	 {
         T    at = a - rhs.GetReal();    // exception guard
         T    bt = b - rhs.GetImag();    // exception guard

		 a = at;
		 b = bt;

		 return(*this);
	 }

	 SIMD_INLINE IOctonion<T>  operator - (const IComplex<T>& rhs) const
	 {
         T    at = a - rhs.GetReal();    // exception guard
         T    bt = b - rhs.GetImag();    // exception guard

		 return IOctonion<T>( at , bt , c , d ,e ,f ,g ,h);
	  }

	 //--------------------------------------------------//

	 SIMD_INLINE IOctonion<T>&  operator -= (const IQuaternion<T>& rhs)
	 {
		 T    at = a - rhs.R_component_1();    // exception guard
		 T    bt = b - rhs.R_component_2();    // exception guard
		 T    ct = c - rhs.R_component_3();    // exception guard
		 T    dt = d - rhs.R_component_4();    // exception guard

		 a = at;
		 b = bt;
		 c = ct;
		 d = dt;

		 return(*this);
	 }

	 SIMD_INLINE IOctonion<T>  operator - (const IQuaternion<T>& rhs) const
	 {
		 T    at = a - rhs.R_component_1();    // exception guard
		 T    bt = b - rhs.R_component_2();    // exception guard
		 T    ct = c - rhs.R_component_3();    // exception guard
		 T    dt = d - rhs.R_component_4();    // exception guard

		 return IOctonion<T>( at , bt , ct , dt ,e ,f ,g ,h);
	 }

	 //--------------------------------------------------//

	 SIMD_INLINE IOctonion<T>&   operator -= (const IOctonion<T>& rhs)
	 {
		 T    at = a - (rhs.R_component_1());    // exception guard
		 T    bt = b - (rhs.R_component_2());    // exception guard
		 T    ct = c - (rhs.R_component_3());    // exception guard
		 T    dt = d - (rhs.R_component_4());    // exception guard
		 T    et = e - (rhs.R_component_5());    // exception guard
		 T    ft = f - (rhs.R_component_6());    // exception guard
		 T    gt = g - (rhs.R_component_7());    // exception guard
		 T    ht = h - (rhs.R_component_8());    // exception guard

		 a = at;
		 b = bt;
		 c = ct;
		 d = dt;
		 e = et;
		 f = ft;
		 g = gt;
		 h = ht;

		 return(*this);
	 }


	 SIMD_INLINE IOctonion<T>   operator - (const IOctonion<T>& rhs) const
	 {
		 T    at = a - (rhs.R_component_1());    // exception guard
		 T    bt = b - (rhs.R_component_2());    // exception guard
		 T    ct = c - (rhs.R_component_3());    // exception guard
		 T    dt = d - (rhs.R_component_4());    // exception guard
		 T    et = e - (rhs.R_component_5());    // exception guard
		 T    ft = f - (rhs.R_component_6());    // exception guard
		 T    gt = g - (rhs.R_component_7());    // exception guard
		 T    ht = h - (rhs.R_component_8());    // exception guard

		 return IOctonion<T>( at , bt , ct , dt ,et ,ft ,gt ,ht);
	 }
	 //--------------------------------------------------//

	 SIMD_INLINE IOctonion<T>&  operator *= (const T& rhs)
	 {
		 T    at = a * rhs;    // exception guard
		 T    bt = b * rhs;    // exception guard
		 T    ct = c * rhs;    // exception guard
		 T    dt = d * rhs;    // exception guard
		 T    et = e * rhs;    // exception guard
		 T    ft = f * rhs;    // exception guard
		 T    gt = g * rhs;    // exception guard
		 T    ht = h * rhs;    // exception guard

		 a = at;
		 b = bt;
		 c = ct;
		 d = dt;
		 e = et;
		 f = ft;
		 g = gt;
		 h = ht;

		 return(*this);
	 }

	 SIMD_INLINE IOctonion<T>  operator * (const T& rhs) const
	 {
		 T    at = a * rhs;    // exception guard
		 T    bt = b * rhs;    // exception guard
		 T    ct = c * rhs;    // exception guard
		 T    dt = d * rhs;    // exception guard
		 T    et = e * rhs;    // exception guard
		 T    ft = f * rhs;    // exception guard
		 T    gt = g * rhs;    // exception guard
		 T    ht = h * rhs;    // exception guard

		 return IOctonion<T>( at , bt , ct , dt ,et ,ft ,gt ,ht);
	 }

	 //--------------------------------------------------//

	 SIMD_INLINE IOctonion<T>&  operator *= (const IComplex<T>& rhs)
	 {
         T    ar = rhs.GetReal();
         T    br = rhs.GetImag();

		 T    at = +a*ar-b*br;
		 T    bt = +a*br+b*ar;
		 T    ct = +c*ar+d*br;
		 T    dt = -c*br+d*ar;
		 T    et = +e*ar+f*br;
		 T    ft = -e*br+f*ar;
		 T    gt = +g*ar-h*br;
		 T    ht = +g*br+h*ar;

		 a = at;
		 b = bt;
		 c = ct;
		 d = dt;
		 e = et;
		 f = ft;
		 g = gt;
		 h = ht;

		 return(*this);
	 }

	 SIMD_INLINE IOctonion<T>  operator * (const IComplex<T>& rhs) const
	 {
         T    ar = rhs.GetReal();
         T    br = rhs.GetImag();

		 T    at = +a*ar-b*br;
		 T    bt = +a*br+b*ar;
		 T    ct = +c*ar+d*br;
		 T    dt = -c*br+d*ar;
		 T    et = +e*ar+f*br;
		 T    ft = -e*br+f*ar;
		 T    gt = +g*ar-h*br;
		 T    ht = +g*br+h*ar;

		 return IOctonion<T>( at , bt , ct , dt ,et ,ft ,gt ,ht);
	 }

	 //--------------------------------------------------//

	 SIMD_INLINE IOctonion<T>&   operator *= (const IQuaternion<T>& rhs)
	 {
		 T    ar = rhs.R_component_1();
		 T    br = rhs.R_component_2();
		 T    cr = rhs.R_component_2();
		 T    dr = rhs.R_component_2();

		 T    at = +a*ar-b*br-c*cr-d*dr;
		 T    bt = +a*br+b*ar+c*dr-d*cr;
		 T    ct = +a*cr-b*dr+c*ar+d*br;
		 T    dt = +a*dr+b*cr-c*br+d*ar;
		 T    et = +e*ar+f*br+g*cr+h*dr;
		 T    ft = -e*br+f*ar-g*dr+h*cr;
		 T    gt = -e*cr+f*dr+g*ar-h*br;
		 T    ht = -e*dr-f*cr+g*br+h*ar;

		 a = at;
		 b = bt;
		 c = ct;
		 d = dt;
		 e = et;
		 f = ft;
		 g = gt;
		 h = ht;

		 return(*this);
	 }

	 SIMD_INLINE IOctonion<T>   operator * (const IQuaternion<T>& rhs) const
	 {
		 T    ar = rhs.R_component_1();
		 T    br = rhs.R_component_2();
		 T    cr = rhs.R_component_2();
		 T    dr = rhs.R_component_2();

		 T    at = +a*ar-b*br-c*cr-d*dr;
		 T    bt = +a*br+b*ar+c*dr-d*cr;
		 T    ct = +a*cr-b*dr+c*ar+d*br;
		 T    dt = +a*dr+b*cr-c*br+d*ar;
		 T    et = +e*ar+f*br+g*cr+h*dr;
		 T    ft = -e*br+f*ar-g*dr+h*cr;
		 T    gt = -e*cr+f*dr+g*ar-h*br;
		 T    ht = -e*dr-f*cr+g*br+h*ar;

		 return IOctonion<T>( at , bt , ct , dt ,et ,ft ,gt ,ht);
	}


	 //--------------------------------------------------//

	 SIMD_INLINE IOctonion<T>&  operator *= (const IOctonion<T>& rhs)
	 {
		 T    ar = (rhs.R_component_1());
		 T    br = (rhs.R_component_2());
		 T    cr = (rhs.R_component_3());
		 T    dr = (rhs.R_component_4());
		 T    er = (rhs.R_component_5());
		 T    fr = (rhs.R_component_6());
		 T    gr = (rhs.R_component_7());
		 T    hr = (rhs.R_component_8());

		 T    at = +a*ar-b*br-c*cr-d*dr-e*er-f*fr-g*gr-h*hr;
		 T    bt = +a*br+b*ar+c*dr-d*cr+e*fr-f*er-g*hr+h*gr;
		 T    ct = +a*cr-b*dr+c*ar+d*br+e*gr+f*hr-g*er-h*fr;
		 T    dt = +a*dr+b*cr-c*br+d*ar+e*hr-f*gr+g*fr-h*er;
		 T    et = +a*er-b*fr-c*gr-d*hr+e*ar+f*br+g*cr+h*dr;
		 T    ft = +a*fr+b*er-c*hr+d*gr-e*br+f*ar-g*dr+h*cr;
		 T    gt = +a*gr+b*hr+c*er-d*fr-e*cr+f*dr+g*ar-h*br;
		 T    ht = +a*hr-b*gr+c*fr+d*er-e*dr-f*cr+g*br+h*ar;

		 a = at;
		 b = bt;
		 c = ct;
		 d = dt;
		 e = et;
		 f = ft;
		 g = gt;
		 h = ht;

		 return(*this);
	 }


	 SIMD_INLINE IOctonion<T>  operator * (const IOctonion<T>& rhs) const
	 {
		 T    ar = (rhs.R_component_1());
		 T    br = (rhs.R_component_2());
		 T    cr = (rhs.R_component_3());
		 T    dr = (rhs.R_component_4());
		 T    er = (rhs.R_component_5());
		 T    fr = (rhs.R_component_6());
		 T    gr = (rhs.R_component_7());
		 T    hr = (rhs.R_component_8());

		 T    at = +a*ar-b*br-c*cr-d*dr-e*er-f*fr-g*gr-h*hr;
		 T    bt = +a*br+b*ar+c*dr-d*cr+e*fr-f*er-g*hr+h*gr;
		 T    ct = +a*cr-b*dr+c*ar+d*br+e*gr+f*hr-g*er-h*fr;
		 T    dt = +a*dr+b*cr-c*br+d*ar+e*hr-f*gr+g*fr-h*er;
		 T    et = +a*er-b*fr-c*gr-d*hr+e*ar+f*br+g*cr+h*dr;
		 T    ft = +a*fr+b*er-c*hr+d*gr-e*br+f*ar-g*dr+h*cr;
		 T    gt = +a*gr+b*hr+c*er-d*fr-e*cr+f*dr+g*ar-h*br;
		 T    ht = +a*hr-b*gr+c*fr+d*er-e*dr-f*cr+g*br+h*ar;

		 return IOctonion<T>( at , bt , ct , dt ,et ,ft ,gt ,ht);
	 }

	 //--------------------------------------------------//

	 SIMD_INLINE IOctonion<T>& operator /= (const T& rhs)
	 {
		 T    at = a / rhs;    // exception guard
		 T    bt = b / rhs;    // exception guard
		 T    ct = c / rhs;    // exception guard
		 T    dt = d / rhs;    // exception guard
		 T    et = e / rhs;    // exception guard
		 T    ft = f / rhs;    // exception guard
		 T    gt = g / rhs;    // exception guard
		 T    ht = h / rhs;    // exception guard

		 a = at;
		 b = bt;
		 c = ct;
		 d = dt;
		 e = et;
		 f = ft;
		 g = gt;
		 h = ht;

		 return(*this);
	 }


	 SIMD_INLINE IOctonion<T> operator / (const T& rhs) const
	 {
		 T    at = a / rhs;    // exception guard
		 T    bt = b / rhs;    // exception guard
		 T    ct = c / rhs;    // exception guard
		 T    dt = d / rhs;    // exception guard
		 T    et = e / rhs;    // exception guard
		 T    ft = f / rhs;    // exception guard
		 T    gt = g / rhs;    // exception guard
		 T    ht = h / rhs;    // exception guard

		 return IOctonion<T>( at , bt , ct , dt ,et ,ft ,gt ,ht);
	}


	 //--------------------------------------------------//

	 SIMD_INLINE IOctonion<T>&  operator /= (const IComplex<T>& rhs)
	 {
         T    ar = rhs.GetReal();
         T    br = rhs.GetImag();

		 T    denominator = ar*ar+br*br;

		 T    at = (+a*ar-b*br)/denominator;
		 T    bt = (-a*br+b*ar)/denominator;
		 T    ct = (+c*ar-d*br)/denominator;
		 T    dt = (+c*br+d*ar)/denominator;
		 T    et = (+e*ar-f*br)/denominator;
		 T    ft = (+e*br+f*ar)/denominator;
		 T    gt = (+g*ar+h*br)/denominator;
		 T    ht = (+g*br+h*ar)/denominator;

		 a = at;
		 b = bt;
		 c = ct;
		 d = dt;
		 e = et;
		 f = ft;
		 g = gt;
		 h = ht;

		 return(*this);
	 }


	 SIMD_INLINE IOctonion<T>  operator / (const IComplex<T>& rhs) const
	 {
         T    ar = rhs.GetReal();
         T    br = rhs.GetImag();

		 T    denominator = ar*ar+br*br;

		 T    at = (+a*ar-b*br)/denominator;
		 T    bt = (-a*br+b*ar)/denominator;
		 T    ct = (+c*ar-d*br)/denominator;
		 T    dt = (+c*br+d*ar)/denominator;
		 T    et = (+e*ar-f*br)/denominator;
		 T    ft = (+e*br+f*ar)/denominator;
		 T    gt = (+g*ar+h*br)/denominator;
		 T    ht = (+g*br+h*ar)/denominator;

		 return IOctonion<T>( at , bt , ct , dt ,et ,ft ,gt ,ht);
	}

	 //--------------------------------------------------//

	 SIMD_INLINE IOctonion<T>&  operator /= (const IQuaternion<T>& rhs)
	 {
		 T    ar = rhs.R_component_1();
		 T    br = rhs.R_component_2();
		 T    cr = rhs.R_component_2();
		 T    dr = rhs.R_component_2();

		 T    denominator = ar*ar+br*br+cr*cr+dr*dr;

		 T    at = (+a*ar+b*br+c*cr+d*dr)/denominator;
		 T    bt = (-a*br+b*ar-c*dr+d*cr)/denominator;
		 T    ct = (-a*cr+b*dr+c*ar-d*br)/denominator;
		 T    dt = (-a*dr-b*cr+c*br+d*ar)/denominator;
		 T    et = (+e*ar-f*br-g*cr-h*dr)/denominator;
		 T    ft = (+e*br+f*ar+g*dr-h*cr)/denominator;
		 T    gt = (+e*cr-f*dr+g*ar+h*br)/denominator;
		 T    ht = (+e*dr+f*cr-g*br+h*ar)/denominator;

		 a = at;
		 b = bt;
		 c = ct;
		 d = dt;
		 e = et;
		 f = ft;
		 g = gt;
		 h = ht;

		 return(*this);
	 }

	 SIMD_INLINE IOctonion<T>  operator / (const IQuaternion<T>& rhs) const
	 {
		 T    ar = rhs.R_component_1();
		 T    br = rhs.R_component_2();
		 T    cr = rhs.R_component_2();
		 T    dr = rhs.R_component_2();

		 T    denominator = ar*ar+br*br+cr*cr+dr*dr;

		 T    at = (+a*ar+b*br+c*cr+d*dr)/denominator;
		 T    bt = (-a*br+b*ar-c*dr+d*cr)/denominator;
		 T    ct = (-a*cr+b*dr+c*ar-d*br)/denominator;
		 T    dt = (-a*dr-b*cr+c*br+d*ar)/denominator;
		 T    et = (+e*ar-f*br-g*cr-h*dr)/denominator;
		 T    ft = (+e*br+f*ar+g*dr-h*cr)/denominator;
		 T    gt = (+e*cr-f*dr+g*ar+h*br)/denominator;
		 T    ht = (+e*dr+f*cr-g*br+h*ar)/denominator;

		 return IOctonion<T>( at , bt , ct , dt ,et ,ft ,gt ,ht);
	 }


	 //--------------------------------------------------//

	 SIMD_INLINE IOctonion<T>&  operator /= (const IOctonion<T>& rhs)
	 {
		 T    ar = (rhs.R_component_1());
		 T    br = (rhs.R_component_2());
		 T    cr = (rhs.R_component_3());
		 T    dr = (rhs.R_component_4());
		 T    er = (rhs.R_component_5());
		 T    fr = (rhs.R_component_6());
		 T    gr = (rhs.R_component_7());
		 T    hr = (rhs.R_component_8());

		 T    denominator = ar*ar+br*br+cr*cr+dr*dr+er*er+fr*fr+gr*gr+hr*hr;

		 T    at = (+a*ar+b*br+c*cr+d*dr+e*er+f*fr+g*gr+h*hr)/denominator;
		 T    bt = (-a*br+b*ar-c*dr+d*cr-e*fr+f*er+g*hr-h*gr)/denominator;
		 T    ct = (-a*cr+b*dr+c*ar-d*br-e*gr-f*hr+g*er+h*fr)/denominator;
		 T    dt = (-a*dr-b*cr+c*br+d*ar-e*hr+f*gr-g*fr+h*er)/denominator;
		 T    et = (-a*er+b*fr+c*gr+d*hr+e*ar-f*br-g*cr-h*dr)/denominator;
		 T    ft = (-a*fr-b*er+c*hr-d*gr+e*br+f*ar+g*dr-h*cr)/denominator;
		 T    gt = (-a*gr-b*hr-c*er+d*fr+e*cr-f*dr+g*ar+h*br)/denominator;
		 T    ht = (-a*hr+b*gr-c*fr-d*er+e*dr+f*cr-g*br+h*ar)/denominator;

		 a = at;
		 b = bt;
		 c = ct;
		 d = dt;
		 e = et;
		 f = ft;
		 g = gt;
		 h = ht;

		 return(*this);
	 }


	 SIMD_INLINE IOctonion<T>  operator / (const IOctonion<T>& rhs) const
	 {
		 T    ar = (rhs.R_component_1());
		 T    br = (rhs.R_component_2());
		 T    cr = (rhs.R_component_3());
		 T    dr = (rhs.R_component_4());
		 T    er = (rhs.R_component_5());
		 T    fr = (rhs.R_component_6());
		 T    gr = (rhs.R_component_7());
		 T    hr = (rhs.R_component_8());

		 T    denominator = ar*ar+br*br+cr*cr+dr*dr+er*er+fr*fr+gr*gr+hr*hr;

		 T    at = (+a*ar+b*br+c*cr+d*dr+e*er+f*fr+g*gr+h*hr)/denominator;
		 T    bt = (-a*br+b*ar-c*dr+d*cr-e*fr+f*er+g*hr-h*gr)/denominator;
		 T    ct = (-a*cr+b*dr+c*ar-d*br-e*gr-f*hr+g*er+h*fr)/denominator;
		 T    dt = (-a*dr-b*cr+c*br+d*ar-e*hr+f*gr-g*fr+h*er)/denominator;
		 T    et = (-a*er+b*fr+c*gr+d*hr+e*ar-f*br-g*cr-h*dr)/denominator;
		 T    ft = (-a*fr-b*er+c*hr-d*gr+e*br+f*ar+g*dr-h*cr)/denominator;
		 T    gt = (-a*gr-b*hr-c*er+d*fr+e*cr-f*dr+g*ar+h*br)/denominator;
		 T    ht = (-a*hr+b*gr-c*fr-d*er+e*dr+f*cr-g*br+h*ar)/denominator;

		 return IOctonion<T>( at , bt , ct , dt ,et ,ft ,gt ,ht);
	 }

	 //----------------------------------------------------------//

	 friend SIMD_INLINE IOctonion<T>  operator - (IOctonion<T> const& o)
	 {
	     return(IOctonion<T>(-o.R_component_1(),
	    		              -o.R_component_2(),
							  -o.R_component_3(),
							  -o.R_component_4(),
							  -o.R_component_5(),
							  -o.R_component_6(),
							  -o.R_component_7(),
							  -o.R_component_8()));
	 }


	 //=========================   Methods =====================================//

     SIMD_INLINE T LengthSquare() const
	  {
	      return a*a + b*b + c*c + d*d + e*e + f*f + g*g + h*h;
	  }

     SIMD_INLINE T Length() const
	  {
          return ISqrt(LengthSquare());
	  }

	  /**
	   * Conjugate IQuaternion
	   */
     SIMD_INLINE IOctonion<T> GetConjugate() const
	  {
		  return   IOctonion<T>(  R_component_1(),
								  -R_component_2(),
								  -R_component_3(),
								  -R_component_4(),
								  -R_component_5(),
								  -R_component_6(),
								  -R_component_7(),
								  -R_component_8());
	  }

	  // Note:    This is the Cayley norm, not the Euclidian norm...
      SIMD_INLINE T Norm() const
	  {
          return ((*this)*(*this).GetConjugate()).GetReal();
	  }

      SIMD_INLINE T Abs() const
	  {
          return(ISqrt( (*this).Norm() ));
	  }



     static SIMD_INLINE IOctonion<T> Exp(IOctonion<T> const & o)
	 {
         T    u = IExp(o.GetReal());
         T    z = (o.GetUnreal().Abs());
		 T    w = ISinc_pi(z);

		 return (IOctonion<T>(ICos(z),
				               w*o.R_component_2(), w*o.R_component_3(),
				               w*o.R_component_4(), w*o.R_component_5(),
				               w*o.R_component_6(), w*o.R_component_7(),
				               w*o.R_component_8()) * u);
	 }


     static SIMD_INLINE IOctonion<T> Cos(IOctonion<T> const & o)
	 {

         T    z =  (o.GetUnreal().Abs());
         T    w = -ISin(o.GetReal())*ISinhc_pi(z);

         return(IOctonion<T>(ICos(o.GetReal())*ICosh(z),
				              w*o.R_component_2(), w*o.R_component_3(),
				              w*o.R_component_4(), w*o.R_component_5(),
				              w*o.R_component_6(), w*o.R_component_7(),
				              w*o.R_component_8()));
	 }


     static SIMD_INLINE IOctonion<T> Sin(IOctonion<T> const & o)
	 {
         T    z = (o.GetUnreal().Abs());
         T    w = +ICos(o.GetReal())*ISinhc_pi(z);

         return(IOctonion<T>(ISin(o.GetReal())*ICosh(z),
	    		              w*o.R_component_2(), w*o.R_component_3(),
				              w*o.R_component_4(), w*o.R_component_5(),
				              w*o.R_component_6(), w*o.R_component_7(),
				              w*o.R_component_8()));
	 }


     static SIMD_INLINE IOctonion<T> Tan(IOctonion<T> const & o)
	 {
         return(Sin(o)/Cos(o));
	 }


     static SIMD_INLINE IOctonion<T> Cosh(IOctonion<T> const & o)
	 {
         return  ((Exp(o)+Exp(-o)) / static_cast<T>(2));
	 }


     static SIMD_INLINE IOctonion<T> Sinh(IOctonion<T> const & o)
	 {
         return((Exp(o)-Exp(-o)) / static_cast<T>(2));
	 }


     static SIMD_INLINE IOctonion<T> Tanh(IOctonion<T>  const & o)
	 {
         return(Sinh(o)/Cosh(o));
	 }


     static SIMD_INLINE IOctonion<T> Pow(IOctonion<T>  const & o, int n)
	 {
		 if(n > 1)
		 {
			 int    m = n>>1;
             IOctonion<T>    result = Pow(o, m);
			 result *= result;
			 if    (n != (m<<1))
			 {
				 result *= o; // n odd
			 }
			 return(result);
		 }
		 else if(n == 1)
		 {
			 return(o);
		 }
		 else if(n == 0)
		 {
			 return(IOctonion<T>(static_cast<T>(1)));
		 }
		 else    /* n < 0 */
		 {
             return(Pow(IOctonion<T>(static_cast<T>(1))/o,-n));
		 }
	}



	 // Note:    There is little point, for the octonions, to introduce the equivalents
	 //            to the complex "arg" and the quaternionic "cylindropolar".
     static SIMD_INLINE IOctonion<T>   Spherical( T const & rho,
												  T const & theta,
												  T const & phi1,
												  T const & phi2,
												  T const & phi3,
												  T const & phi4,
												  T const & phi5,
												  T const & phi6)
	 {

		 //T    a = cos(theta)*cos(phi1)*cos(phi2)*cos(phi3)*cos(phi4)*cos(phi5)*cos(phi6);
		 //T    b = sin(theta)*cos(phi1)*cos(phi2)*cos(phi3)*cos(phi4)*cos(phi5)*cos(phi6);
		 //T    c = sin(phi1)*cos(phi2)*cos(phi3)*cos(phi4)*cos(phi5)*cos(phi6);
		 //T    d = sin(phi2)*cos(phi3)*cos(phi4)*cos(phi5)*cos(phi6);
		 //T    e = sin(phi3)*cos(phi4)*cos(phi5)*cos(phi6);
		 //T    f = sin(phi4)*cos(phi5)*cos(phi6);
		 //T    g = sin(phi5)*cos(phi6);
		 //T    h = sin(phi6);

		 T    courrant = static_cast<T>(1);
         T    h = ISin(phi6);courrant *= ICos(phi6);
         T    g = ISin(phi5)*courrant; courrant *= ICos(phi5);
         T    f = ISin(phi4)*courrant; courrant *= ICos(phi4);
         T    e = ISin(phi3)*courrant; courrant *= ICos(phi3);
         T    d = ISin(phi2)*courrant; courrant *= ICos(phi2);
         T    c = ISin(phi1)*courrant; courrant *= ICos(phi1);

         T    b = ISin(theta)*courrant;
         T    a = ICos(theta)*courrant;

		 return(rho*IOctonion<T>(a,b,c,d,e,f,g,h));
	 }


     static SIMD_INLINE  IOctonion<T> Multipolar( T const & rho1,
												  T const & theta1,
												  T const & rho2,
												  T const & theta2,
												  T const & rho3,
												  T const & theta3,
												  T const & rho4,
												  T const & theta4)
	 {

         T    a = rho1*ICos(theta1);
         T    b = rho1*ISin(theta1);
         T    c = rho2*ICos(theta2);
         T    d = rho2*ISin(theta2);
         T    e = rho3*ICos(theta3);
         T    f = rho3*ISin(theta3);
         T    g = rho4*ICos(theta4);
         T    h = rho4*ISin(theta4);

		 return(IOctonion<T>(a,b,c,d,e,f,g,h));
	 }


     static SIMD_INLINE IOctonion<T> Cylindrical( T const & r,
												  T const & angle,
												  T const & h1,
												  T const & h2,
												  T const & h3,
												  T const & h4,
												  T const & h5,
												  T const & h6)
	 {

         T    a = r*ICos(angle);
         T    b = r*ISin(angle);

		 return(IOctonion<T>(a,b,h1,h2,h3,h4,h5,h6));
	 }


	 //----------[ output operator ]----------------------------
	  /**
	  * Provides output to standard output stream.
	  */
	  friend std::ostream& operator <<(std::ostream& oss, const IOctonion<T>& q)
	  {
		  oss << "(" << "Re: " << q.a << " Im: " << "[" << q.b << "," << q.c << "," << q.d << "," << q.e << ","
				                                        << q.f << "," << q.g << "," << q.h << "] )";
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


	  /**
	      * The multiplicitive identity Octonion
	      */
	     static const IOctonion<T> IDENTITY;
	     /**
	      * The additive identity Octonion.
	      */
	     static const IOctonion<T> ZERO;

};

template<class T> const IOctonion<T> IOctonion<T>::IDENTITY(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
template<class T> const IOctonion<T> IOctonion<T>::ZERO(0.0, 0.0, 0.0, 0.0, 0.0, 0.0 , 0.0, 0.0);

//--------------------------------------
// Typedef shortcuts for Octonion
//-------------------------------------

using IOctonionr    = IOctonion<Real>;
using IOctonionf    = IOctonion<float>;
using IOctoniond    = IOctonion<double>;
using IOctonioni    = IOctonion<std::int32_t>;
using IOctonionui   = IOctonion<std::uint32_t>;
using IOctonionb    = IOctonion<std::int8_t>;
using IOctonionub   = IOctonion<std::uint8_t>;


} /* namespace */


#endif /* IOCTONION_HPP_ */
