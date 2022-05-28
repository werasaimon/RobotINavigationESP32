 /********************************************************************************
 *
 * IComplex.h
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
///

#ifndef ICOMPLEX_H_
#define ICOMPLEX_H_


#include "IMatrix2x2.h"


namespace IMath
{

// define Complex class
//=========================================================================
template<class T> class  IComplex
{
 private:

	T real,imag;


 public:

	// default constructor
	SIMD_INLINE IComplex(void)
    : real(0.0) ,
      imag(0.0)
    {}

	// constructor
	SIMD_INLINE IComplex(T r)
	: real(r) , imag(0.0)
	{}

	// constructor
	SIMD_INLINE IComplex(T r,T i)
	: real(r) , imag(i)
	{}



	// return the real part
    SIMD_INLINE T GetReal(void) const { return(this->real); }

	// return the imagimary part
    SIMD_INLINE T GetImag(void) const { return(this->imag); }


	// return the phase argument
    SIMD_INLINE T GetArgument(void) const { return atan2(this->imag,this->real); }

	// return the magnitude square
    SIMD_INLINE T LengthSquare(void) const { return (real*real+imag*imag); }

	// return the magnitude
    SIMD_INLINE T Length(void) const { return ISqrt(real*real+imag*imag); }


	// return the real matrix
    SIMD_INLINE IMatrix2x2<T> GetMatrix() const
	{
	  return IMatrix2x2<T> ( real , -imag ,
			                 imag ,  real );
	}


	// return the complex conjugate
    SIMD_INLINE IComplex<T> GetConjugate(void) { return IComplex<T>(this->real,-this->imag); }


	// overload the + operator to add 2 complex numbers
	SIMD_INLINE IComplex<T> operator+(IComplex<T> z) const
	{
		return IComplex<T>(this->real+z.real,this->imag+z.imag);
	}

	// overload the + operator to add complex and a T
	SIMD_INLINE IComplex<T> operator+(T a) const
	{
		return IComplex<T>(this->real+a,this->imag);
	}

	// This is a friend function
	// overload the + operator to add a T and a complex
	friend SIMD_INLINE IComplex<T> operator+(T a, IComplex<T> z)
	{
		return IComplex<T>(a+z.real,z.imag);
	}

	// overload the - operator to subtract 2 complex numbers
	SIMD_INLINE IComplex<T> operator-(IComplex<T> z) const
	{
		return IComplex<T>(this->real-z.real,this->imag-z.imag);
	}

	// overload the - operator to subtract a T from a complex
	SIMD_INLINE IComplex<T> operator-(T a) const
	{
		return IComplex<T>(this->real-a,this->imag);
	}

	// overload the - operator to subtract a complex from a T
	friend SIMD_INLINE IComplex<T> operator-(T a, IComplex<T> z)
	{
		return IComplex<T>(a-z.real,-z.imag);
	}

	// overload the - operator to take the negative of a complex
	friend SIMD_INLINE IComplex<T> operator-(IComplex<T> z)
	{
		return IComplex<T>(-z.real,-z.imag);
	}

	// overload the * operator to multiply two complex numbers
	SIMD_INLINE IComplex<T> operator*(IComplex<T> z) const
	{
		return IComplex<T>(this->real*z.real-this->imag*z.imag,
				            this->real*z.imag+this->imag*z.real);
	}

	// overload the * operator to multiply a complex by a T
	SIMD_INLINE IComplex<T> operator*(T a) const
	{
		return IComplex<T>(real*a,imag*a);
	}

	// overload the * operator to multiply a T by a complex
	friend SIMD_INLINE IComplex<T> operator*(T a, IComplex<T> z)
	{
		return IComplex<T>(a*z.real,a*z.imag);
	}

	// overload the / operator to divide two complex numbers
	SIMD_INLINE IComplex<T> operator/(IComplex<T> z) const
	{
        IComplex<T> top((*this)*z.GetConjugate());
        T bottom(z.LengthSquare());
		IComplex<T> quo(top/bottom);
		return quo;
	}

	// overload the / operator to divide a complex number by a T
	SIMD_INLINE IComplex<T> operator/(T a) const
	{
		return IComplex<T>(this->real/a,this->imag/a);
	}

	// overload the / operator to divide a T by a complex
	friend SIMD_INLINE IComplex<T> operator/(T a, IComplex<T> z)
	{
        IComplex<T> top((a)*z.GetConjugate());
        T bottom(z.LengthSquare());
		IComplex<T> quo(top/bottom);
		return quo;
	}

	// overload the += operator
	SIMD_INLINE const IComplex<T>& operator+=(const IComplex<T>& z)
    {
		this->real+=z.real;
		this->imag+=z.imag;
		return *this;
	}

	// overload the -= operator
	SIMD_INLINE const IComplex<T>& operator-=(const IComplex<T>& z)
	{
		this->real-=z.real;
		this->imag-=z.imag;
		return *this;
	}

	// overload the == operator
	SIMD_INLINE bool operator==(IComplex<T> z)
	{
		if (this->real == z.real &&
            this->imag == z.imag)
		{
			return true;
		}
		else
		{
			return false;
		}
	}


    //--------------------------//

	static const IComplex<T> j;
	static const IComplex<T> i;

	//--------------------------//

	// take the square root of a complex number
    static SIMD_INLINE IComplex<T> Sqrt(IComplex<T> z)
    {
		T zsqre,zsqim;

        zsqre = ISqrt(0.5*(z.Length()+z.GetReal()));
        zsqim = ISqrt(0.5*(z.Length()-z.GetImag()));

        if (z.GetImag() >= 0.0)
		{
			return IComplex<T>(zsqre,zsqim);
		}
		else
		{
			return IComplex<T>(zsqre,-zsqim);
		}
	}

	// take the natural log of a complex number
    static SIMD_INLINE IComplex<T> Log(IComplex<T> z)
	{
		if (z.real < 0 && z.imag == 0.0)
		{
            return ILog(z.Length())+j*M_PI;
		}
		else
		{
            return ILog(z.Length())+j*z.GetArgument();
		}
	}

	//========= Triqonometriya Function ==========//

	// raise e to a complex number
    static SIMD_INLINE IComplex<T> Exp(IComplex<T> z)
	{
		return IExp(z.real)*(ICos(z.imag)+j*ISin(z.imag));
	}

	// raise a complex number to a T
    static SIMD_INLINE IComplex<T> Pow(IComplex<T> z, T c)
	{
        return Exp(c*log(z));
	}


	// take the sin of a complex number
    static SIMD_INLINE IComplex<T> Sin(IComplex<T> z)
	{
        return 0.5*(-j)*Exp(j*z)-0.5*(-j)*Exp(-j*z);
	}

	// take the cos of a complex number
    static SIMD_INLINE IComplex<T> Cos(IComplex<T> z)
	{
        return 0.5*Exp(j*z)+0.5*Exp(-j*z);
	}

	// take the tan of a complex number
    static SIMD_INLINE IComplex<T> Tan(IComplex<T> z) { return Sin(z)/Cos(z); }

	// take the sec of a complex number
    static SIMD_INLINE IComplex<T> Sec(IComplex<T> z) { return 1/Cos(z); }

	// take the csc of a complex number
    static SIMD_INLINE IComplex<T> Csc(IComplex<T> z) { return 1/Sin(z); }

	// take the cot of a complex number
    static SIMD_INLINE IComplex<T> Cot(IComplex<T> z) { return Cos(z)/Sin(z); }

	// take the sinh of a complex number
    static SIMD_INLINE IComplex<T> Sinh(IComplex<T> z) { return (Exp(z)-Exp(-z))/2.0; }

	// take the cosh of a complex number
    static SIMD_INLINE IComplex<T> Cosh(IComplex<T> z) { return (Exp(z)+Exp(-z))/2.0; }

	// take the tanh of a complex number
    static SIMD_INLINE IComplex<T> Tanh(IComplex<T> z) { return Sinh(z)/Cosh(z); }

	// take the sech of a complex number
    static SIMD_INLINE IComplex<T> Sech(IComplex<T> z) { return 1/Cosh(z); }

	// take the csch of a complex number
    static SIMD_INLINE IComplex<T> Csch(IComplex<T> z) { return 1/Sinh(z); }

	// take the coth of a complex number
    static SIMD_INLINE IComplex<T> Coth(IComplex<T> z) { return Cosh(z)/Sinh(z); }

	// take the asin of a complex number
    static SIMD_INLINE IComplex<T> Asin(IComplex<T> z) { return -j*Log(j*z+Sqrt(1.0-z*z)); }

	// take the acos of a complex number
    static SIMD_INLINE IComplex<T> Acos(IComplex<T> z) { return -j*Log(z+Sqrt(z*z-1.0)); }

	// take the atan of a complex number
    static SIMD_INLINE IComplex<T> Atan(IComplex<T> z) { return (0.5*j)*Log((j+z)/(j-z)); }

	// take the asinh of a complex number
    static SIMD_INLINE IComplex<T> Asinh(IComplex<T> z) { return Log(z+Sqrt(z*z+1.0)); }

	// take the acosh of a complex number
    static SIMD_INLINE IComplex<T> Acosh(IComplex<T> z) { return Log(z+Sqrt(z*z-1.0)); }

	// take the atanh of a complex number
    static SIMD_INLINE IComplex<T> Atanh(IComplex<T> z) { return 0.5*Log((1.0+z)/(1.0-z)); }



	// round a complex number
    SIMD_INLINE IComplex<T> Rnd(int precision)
	{
		T rnum,inum;
		int tnum;

		rnum = this->real*pow(10,precision);
        tnum = (rnum < 0 ? rnum-0.5 : rnum + 0.5);
		rnum = tnum/pow(10,precision);

		inum = this->imag*pow(10,precision);
        tnum = (inum < 0 ? inum-0.5 : inum + 0.5);
		inum = tnum/pow(10,precision);

		return IComplex<T>(rnum,inum);
	}

    static SIMD_INLINE IComplex<T> Rnd(IComplex<T> z, int precision)
	{
		T rnum,inum;
		int tnum;

		rnum = z.real*pow(10,precision);
        tnum = (rnum < 0 ? rnum-0.5 : rnum + 0.5);
		rnum = tnum/pow(10,precision);

		inum = z.imag*pow(10,precision);
        tnum = (inum < 0 ? inum-0.5 : inum + 0.5);
		inum = tnum/pow(10,precision);


		return IComplex<T>(rnum,inum);
	}

	// end of Complex class
	//=========================================================================

	// round a number
    static SIMD_INLINE T Rnd(T num, int precision)
	{
		T rnum;
		int tnum;

		rnum = num*pow(10,precision);
        tnum = (rnum < 0 ? rnum-0.5 : rnum + 0.5);
		rnum = tnum/pow(10,precision);

		return rnum;
	}


	// create an inserter function so
	// complex types work with cout
	friend std::ostream& operator<<(std::ostream& stream, IComplex<T> z)
	{
        stream << "(" << "Re: " << z.real << " Im: " << z.imag << ")";

		return stream;
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


template<class T> const IComplex<T> IComplex<T>::i = IComplex<T>(0.0f, 1.0f);
template<class T> const IComplex<T> IComplex<T>::j = IComplex<T>(0.0f, 1.0f);


//--------------------------------------
// Typedef shortcuts for Complex
//-------------------------------------

using IComplexr    = IComplex<Real>;
using IComplexf    = IComplex<float>;
using IComplexd    = IComplex<double>;
using IComplexi    = IComplex<std::int32_t>;
using IComplexui   = IComplex<std::uint32_t>;
using IComplexb    = IComplex<std::int8_t>;
using IComplexub   = IComplex<std::uint8_t>;

} /* namespace */




#endif /* SRC_MATHS_ICOMPLEX_H_ */
