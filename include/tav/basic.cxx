/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file basic.cxx Basic templated function in #tav namespace
*/
#ifndef _BASIC_H_
#define _BASIC_H_
#include <debug.h> //FatalError depends on this
/*
Convention to use capital letter for all first letters in full words
10/17/2005 Basics Read Write to file with numbers
           Swap function call by reference
11/09/2005 CopySign -double template
11/17/2005 Even, Odd, Round
01/07/2006 include fstream, make simple input/output operations in scope 
03/13/2006 New basic stream operators << and >> are given by || which ignores 
         comment lines
	 do not forget "using tav::operator ||;"
	 operator || is simplified using peek
03/14/2006 operator (char << variable)  are introduced append mode 
           < overwrite mode, operaotr only works with derived classesx
	   
	   TestFile(char*) is introduced
04/14/2006 Pow(Template, int) int>1
04/18/2006 disable Read, because it conflicts with Read<tensor>
07/16/2008 In TestFile, must use const
09/09/2009 OppositeSign test of two numbers
*/


/*! \namespase tav
\brief Templated mathematical library of functions and containers
 */
namespace tav{


// even and odd
//! \brief round a templated float number to an integer 
template <class cnumber> 
inline int	Round (cnumber d) { return (d>0) ? (int)(d+0.5) : -(int)(-d+0.5); }


template<class inumber> 
inline inumber Even (inumber x) {return!(x&1);}
template<class inumber> 
inline inumber Odd (inumber x) {return (x&1);}
	

// ADDITIONAL STANDARD FUNCTIONS AS TEMPLATES



// return sign of number


template <class T> inline int Sign (T x)
{ return (x > 0) ? 1 : (x < 0) ? -1 : 0; }
template <class T> inline T Abs (T x)
  { return (x >= 0) ? x : -x; } // may conflict with complex double=Abs(complex)
template <class T> inline int Parity (T x)
{ return (1-(((x)&1)<<1));}
  //this returns parity of the number even +1, odd -1 only integer 

// return first number with sign of second number
//Copy sign returns x having the same sign as y
template <class T, class Q> inline T CopySign (T x, Q y)
{ return (y < 0) ? ((x < 0) ? x : -x) : ((x > 0) ? x : -x); }
//!\breif Check if number have the opposite sign 0-same 1-opposite
template <class T, class Q> inline T OppositeSign (T x, Q y)
{ return (x<0)^(y<0); }
//-----------------------------------------------------------------------------//
// MpMin(), MpMin(): min and max templates for upto 9 arguments
//-----------------------------------------------------------------------------//

template <class T> inline T Min (T x, T y) 
{ return (x<y)?x:y; }
template <class T> inline T Min (T x, T y, T z) 
{ return (x<y)?((x<z)?x:z):((y<z)?y:z); }
template <class T> inline T Min (T x, T y, T z, T w) 
{ return Min(Min(x,y),Min(z,w)); }
template <class T> inline T Min (T x1,T x2,T x3,T x4,T x5) 
{ return Min(Min(x1,x2),Min(x3,x4),x5); }
template <class T> inline T Min (T x1,T x2,T x3,T x4,T x5,T x6) 
{ return Min(Min(x1,x2,x3),Min(x4,x5,x6)); }
template <class T> inline T Min (T x1,T x2,T x3,T x4,T x5,T x6,T x7) 
{ return Min(Min(x1,x2,x3),Min(x4,x5,x6),x7); }
template <class T> inline T Min (T x1,T x2,T x3,T x4,T x5,T x6,T x7,T x8) 
{ return Min(Min(x1,x2,x3),Min(x4,x5,x6),Min(x7,x8)); }
template <class T> inline T Min (T x1,T x2,T x3,T x4,T x5,T x6,T x7,T x8, T x9) 
{ return Min(Min(x1,x2,x3),Min(x4,x5,x6),Min(x7,x8,x9)); }

template <class T> inline T Max (T x, T y) 
{ return (x>y)?x:y; }
template <class T> inline T Max (T x, T y, T z) 
{ return (x>y)?((x>z)?x:z):((y>z)?y:z); }
template <class T> inline T Max (T x, T y, T z, T w) 
{ return Max(Max(x,y),Max(z,w)); }
template <class T> inline T Max (T x1,T x2,T x3,T x4,T x5) 
{ return Max(Max(x1,x2),Max(x3,x4),x5); }
template <class T> inline T Max (T x1,T x2,T x3,T x4,T x5,T x6) 
{ return Max(Max(x1,x2,x3),Max(x4,x5,x6)); }
template <class T> inline T Max (T x1,T x2,T x3,T x4,T x5,T x6,T x7) 
{ return Max(Max(x1,x2,x3),Max(x4,x5,x6),x7); }
template <class T> inline T Max (T x1,T x2,T x3,T x4,T x5,T x6,T x7,T x8) 
{ return Max(Max(x1,x2,x3),Max(x4,x5,x6),Max(x7,x8)); }
template <class T> inline T Max (T x1,T x2,T x3,T x4,T x5,T x6,T x7,T x8, T x9) 
{ return Max(Max(x1,x2,x3),Max(x4,x5,x6),Max(x7,x8,x9)); }

//-----------------------------------------------------------------------------//
// raise to n-th power for n = 2 to 12
//-----------------------------------------------------------------------------//

template <class T> inline T    Sqr  (T x) { return x*x; }
template <class T> inline T   Cube  (T x) { return x*x*x; }
template <class T> inline T Pow2nd  (T x) { return x*x; }
template <class T> inline T Pow3rd  (T x) { return x*x*x; }
template <class T> inline T Pow4th  (T x) { T y = x*x; return y*y; }
template <class T> inline T Pow5th  (T x) { T y = x*x; return y*y*x; }
template <class T> inline T Pow6th  (T x) { T y = x*x*x; return y*y; }
template <class T> inline T Pow7th  (T x) { T y = x*x*x; return y*y*x; }
template <class T> inline T Pow8th  (T x) { T y = x*x; y *= y; return y*y; }
template <class T> inline T Pow9th  (T x) { T y = x*x; y *= y; return y*y*x; }
template <class T> inline T Pow10th (T x) { T y = x*x; T z = y*y; return z*z*y; }
template <class T> inline T Pow11th (T x) { T y = x*x; T z = y*y; return z*z*y*x; }
template <class T> inline T Pow12th (T x) { T y = x*x*x; y *= y; return y*y; }
template <class T> T Pow(T x, unsigned int n) {T y=x; 
for (int i=1;i<n;i++) y*=x; 
  return y;
}
// Direct swap of A and B if you need swap pointers use Swap(&A,&B)
template <class T> inline void Swap (T &a, T &b)
{ T c(a); a = b; b = c; }


}
#endif

