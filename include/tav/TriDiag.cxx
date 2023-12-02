/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file TriDiag.cxx
\brief Diagonalization and linear equations with tri-diagonal matrices
\ingroup gp_nr
 */
 /*
 2/26/09	Templated names, call by template with reference 
			fabs is changed to abs, 
			2.0*x is changed as a sum x+x
2/27/09     Remove Sqr and dependence on basic
 */
#ifndef _TRIDIAG_CXX_
#define _TRIDIAG_CXX_
//#include "basic.cxx" //uses Sqr
//#include "debug.h" //uses FatalError
#include <iostream>
#define FatalError_TriDiag(message)              { std::cerr <<"Error: "<< message <<std::endl<< __FILE__ << ":" << __LINE__ << std::endl; exit(1); }
#include <cmath> //uses abs for double
using std::abs;
namespace tav{
/*!\brief Sove linear equation with tri-diagonal matrix

\f[
\left[\begin{array}{ccccc}
b_{0} & c_{0}\\
a_{1} & b_{1} & c_{1}\\
 & a_{2} & b_{2} & c_{2}\\
 &  & . & . & .\\
 &  &  & a_{n-1} & b_{n-1}\end{array}\right]\left[\begin{array}{c}
u_{0}\\
u_{1}\\
u_{2}\\
.\\
u_{n-1}\end{array}\right]=\left[\begin{array}{c}
r_{0}\\
r_{1}\\
r_{2}\\
.\\
r_{n-1}\end{array}\right]
\f]
*/
template <class dclass> 
void SolveTriDiag(
unsigned int n, ///<dimension of the matrix
dclass &a, ///<left off-diagonal a[0] not used
dclass &b, ///<diagonal 
dclass &c, ///< right off-diagonal c[n-1] not used
dclass &r, ///< right-hand side 
dclass &u ///< output: return solution
) 
//a from 0 b from 0 c from 1 numbering by row
{
	int j;
//	typedef typeof(a[0]) cnumber;
#if __cplusplus <= 199711L
    typedef typeof(a[0]) cnumber;
#else
    auto xxxc=a[0];
    typedef decltype(xxxc) cnumber;
#endif
	cnumber bet;
	//int n=a.size();
	cnumber *gam= new cnumber[n];
	if (b[0] == 0.0) std::cerr<<"Error 1 in tridag\n";
	u[0]=r[0]/(bet=b[0]);
	for (j=1;j<n;j++) {
		gam[j]=c[j-1]/bet;
		bet=b[j]-a[j]*gam[j];
		if (bet == 0.0) std::cerr<<"Error 2 in tridag\n";
		u[j]=(r[j]-a[j]*u[j-1])/bet;
	}
	for (j=(n-2);j>=0;j--)
		u[j] -= gam[j+1]*u[j+1];
	delete [] gam;
};
  //this will work for complex as long as we have all stuff well defined
  //dclass is more general
//#define pythag(a,b) sqrt(a*a+b*b)
/**\brief Pythagorian length \f$ \sqrt(a^2+b^2) \f$
*/
template <class dclass, class qclass> 
dclass pythag(const dclass& a, const qclass& b)
{
  //return sqrt(a*a+b*b); //this is what it does
 
	dclass absa,absb;
	absa=abs(a);
	absb=abs(b);
	if (absa > absb) 
	{dclass frac=absb/absa; return absa*sqrt(1.0+frac*frac);}
	if (absb==0.0) return 0.0; 
	dclass frac=absa/absb; 
	return absb*sqrt(1.0+frac*frac);
	}



#ifndef SIGN
//#define ABS(a) ((a)>=0.0 ? a: -a )
#define SIGN(a,b) ((b) >= 0.0 ? ((a)>=0.0 ? a : -a ) : ((a)>=0.0 ? -a : a )) 
// above will not work for complex
//#define SIGN(a,b) ((b) >= 0.0 ? abs(a) : -abs(a))
//at present b>=0.0 is always true for complex
//#define SIGN(a,b) abs(a)*b/abs(b) //complex require somparison for real part
//only effects shift which will not be very successful if above is bad
//template <class dclass> 
/*! \brief Eigenvalues of tri-diagonal matrix
This routine presently is not compatible with complex, but can work potentially
-The routine is compatible with multipole precision arprec library.   
\return eigenvalues in vector d[]
*/
template <class cvector1, class cvector2>
void Tqli (
int n, ///< dimension of the matrix 
cvector1 &d, ///< input: diagonal values, access d[]; output-eigenvalues
cvector2 &e ///< off-diagonal \f$ v_n=M_{n,n+1} \f$, v[n-1] has no meaning. 
)
  //d diagonal, e off-diagonal, e[n-1] has no meaning but must exist
  //, n size of array
//void tqli(float d[], float e[], int n, float **z)
{
  int m,l,iter,i,k;
//  typedef typeof(e[0]) cnumber; //generally complex
//  typedef typeof(d[0]) dnumber; //generally double
#if __cplusplus <= 199711L
    typedef typeof(e[0]) cnumber;
    typedef typeof(d[0]) dnumber;
#else
    auto xxxc=e[0];
    auto xxxd=d[0];
    typedef decltype(xxxc) cnumber;
    typedef decltype(xxxd) dnumber;
#endif
	cnumber s,r,p,g,f,dd,c,b;

	//	for (i=1;i<n;i++) e[i-1]=e[i];
	e[n-1]=0.0;
	for (l=0;l<n;l++) {
		iter=0;
		do {
			for (m=l;m<n-1;m++) {
				dd=abs(d[m])+abs(d[m+1]);
				if ((abs(e[m])+dd) == dd) break;
			}
			if (m != l) {
			  if (iter++ == 30) FatalError_TriDiag("Too many iterations in tqli");
				g=(d[l+1]-d[l])/(e[l]+e[l]);
				r=pythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=pythag(f,g));
					if (r == 0.0) {
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					//r=(d[i]-g)*s+2.0*c*b;
					r=c*b; r+=r+(d[i]-g)*s;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
}
/*! \brief Eigenvalues and eigenvectors of tri-diagonal matrix 
\return eigenvalues in vector d[] and egenvectors in z
*/
template <class cmatrix, class cvector1, class cvector2>
int Tqli (int n, ///<dimension of the matrix
cvector1 &d, ///< input: diagonal values, access d[]; output-eigenvalues
cvector2 &e, ///< off-diagonal \f$ v_n=M_{n,n+1} \f$, v[n-1] has no meaning.
cmatrix &z ///<input unit matrix, access z[][]; output eigenvalues columns
)
  //d diagonal, e off-diagonal, e[n-1] has no meaning but must exist
  //, n size of array z identity matrix
//void tqli(float d[], float e[], int n, float **z)
{
  int m,l,iter,i,k;
//   typedef typeof(e[0]) cnumber; //generally complex
//  typedef typeof(d[0]) dnumber; //generally double
#if __cplusplus <= 199711L
    typedef typeof(e[0]) cnumber;
    typedef typeof(d[0]) dnumber;
#else
    auto xxxc=e[0];
    auto xxxd=d[0];
    typedef decltype(xxxc) cnumber;
    typedef decltype(xxxd) dnumber;
#endif
	cnumber s,r,p,g,f,dd,c,b;


	//	for (i=1;i<n;i++) e[i-1]=e[i];
	e[n-1]=0.0;
	for (l=0;l<n;l++) {
		iter=0;
		do {
			for (m=l;m<n-1;m++) {
				dd=abs(d[m])+abs(d[m+1]);
				if ((abs(e[m])+dd) == dd) break;
			}
			if (m != l) {
			  if (iter++ == 30) FatalError_TriDiag("Too many iterations in tqli");
			  				g=(d[l+1]-d[l])/(e[l]+e[l]);
				r=pythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=pythag(f,g));
					if (r == 0.0) {
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=c*b; r+=r+(d[i]-g)*s;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
				for (k=0;k<n;k++) {
						f=z[k][i+1];
						z[k][i+1]=s*z[k][i]+c*f;
						z[k][i]=c*z[k][i]-s*f;
					}
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
return 1;
}
/*! \brief Eigenvalues and eigenvectors of tri-diagonal matrix, transposed form 
\return eigenvalues in vector d[] and egenvectors in z (transposed form
*/
template <class cmatrix, class cvector1, class cvector2>
int TqliT (int n, ///< dimension of the matrix
cvector1 &d, ///< input: diagonal values, access d[]; output-eigenvalues
cvector2 &e, ///< off-diagonal \f$ v_n=M_{n,n+1} \f$, v[n-1] has no meaning.
cmatrix &z///< input: unit matrix; output: eigenvalues in rows
)
  //d diagonal, e off-diagonal, e[n-1] has no meaning but must exist
  //, n size of array z identity matrix
//void tqli(float d[], float e[], int n, float **z)
    //return 1 ok and 0 too many iterations
{
  int m,l,iter,i,k;
//  typedef typeof(e[0]) cnumber; //generally complex
//  typedef typeof(d[0]) dnumber; //generally double
#if __cplusplus <= 199711L
    typedef typeof(e[0]) cnumber;
    typedef typeof(d[0]) dnumber;
#else
    auto xxxc=e[0];
    auto xxxd=d[0];
    typedef decltype(xxxc) cnumber;
    typedef decltype(xxxd) dnumber;
#endif

	cnumber s,r,p,g,f,dd,c,b;

	//	for (i=1;i<n;i++) e[i-1]=e[i];
	e[n-1]=0.0;
	for (l=0;l<n;l++) {
		iter=0;
		do {
			for (m=l;m<n-1;m++) {
				dd=abs(d[m])+abs(d[m+1]);
				if ((abs(e[m])+dd) == dd) break;
			}
			if (m != l) {
			  if (iter++ == 30) FatalError_TriDiag("Too many iterations in tqli");
			  				g=(d[l+1]-d[l])/(e[l]+e[l]);
				r=pythag(g,cnumber(1.0));
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=pythag(f,g));
					if (r == 0.0) {
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=c*b; r+=r+(d[i]-g)*s;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
				for (k=0;k<n;k++) {
				  //f=z[k][i+1];
				  f=z[i+1][k];
				  //z[k][i+1]=s*z[k][i]+c*f;
				  z[i+1][k]=s*z[i][k]+c*f;
				  //z[k][i]=c*z[k][i]-s*f;
				  z[i][k]=c*z[i][k]-s*f;
					}
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
return 1;
}


#undef SIGN 
#endif
}
#endif

