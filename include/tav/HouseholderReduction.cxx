/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/** 
\author Alexander Volya  <http://www.volya.net>
\file HouseholderReduction.cxx
\brief Householder Reduction: Reduce general matrix to tri-diagonal
\ingroup gp_nr
*/
#ifndef _HOUSEHOLDERREDUCTION_CXX_
#define _HOUSEHOLDERREDUCTION_CXX_
#include <cmath> //uses abs for double
using std::abs;
namespace tav {
/** \brief Householder reduction 
*/
  template <class cmatrix, class cvector1, class cvector2>
  void HouseholderReduction(int n, ///<dimension 
cmatrix &a, ///< real symmetric matrix, output is replaced by orthogonal matrix \e Q
cvector1 &d, ///< diagonal elements 
cvector2 &e ///< off-diagonal elements \e e[n-1]=0 and is irrelevant
)  
{
	int l,k,j,i;
//	typedef typeof(a[0][0]) cnumber;
#if __cplusplus <= 199711L
    typedef typeof(a[0][0]) cnumber;
#else
    auto xxxc=a[0][0];
    typedef decltype(xxxc) cnumber;
#endif
	cnumber scale,hh,h,g,f;
	//	int n=d.size();
	for (i=n-1;i>0;i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 0) {
			for (k=0;k<l+1;k++)
				scale += abs(a[i][k]);
			if (scale == 0.0)
				e[i]=a[i][l];
			else {
				for (k=0;k<l+1;k++) {
					a[i][k] /= scale;
					h += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				a[i][l]=f-g;
				f=0.0;
				for (j=0;j<l+1;j++) {
				// Next statement can be omitted if eigenvectors not wanted
					a[j][i]=a[i][j]/h;
					g=0.0;
					for (k=0;k<j+1;k++)
						g += a[j][k]*a[i][k];
					for (k=j+1;k<l+1;k++)
						g += a[k][j]*a[i][k];
					e[j]=g/h;
					f += e[j]*a[i][j];
				}
				hh=f/(h+h);
				for (j=0;j<l+1;j++) {
					f=a[i][j];
					e[j]=g=e[j]-hh*f;
					for (k=0;k<j+1;k++)
						a[j][k] -= (f*e[k]+g*a[i][k]);
				}
			}
		} else
			e[i]=a[i][l];
		d[i]=h;
	}
	// Next statement can be omitted if eigenvectors not wanted
	d[0]=0.0;
	e[0]=0.0;
	// Contents of this loop can be omitted if eigenvectors not
	//	wanted except for statement d[i]=a[i][i];
	for (i=0;i<n;i++) {
		l=i;
		if (d[i] != 0.0) {
			for (j=0;j<l;j++) {
				g=0.0;
				for (k=0;k<l;k++)
					g += a[i][k]*a[k][j];
				for (k=0;k<l;k++)
					a[k][j] -= g*a[k][i];
			}
		}
		d[i]=a[i][i];
		a[i][i]=1.0;
		for (j=0;j<l;j++) a[j][i]=a[i][j]=0.0;
	}
	for (i=1;i<n;i++) e[i-1]=e[i];  //create return so offdiagonal 0,1,..
	e[n-1]=0.0;

}
/*! \brief Housholder Method, reduction of a symmetic matrix to a tridiagonal form. Output is in transposed form. 
*/
  template <class cmatrix, class cvector1, class cvector2>
  void HouseholderReductionT(
      int n,///< dimension of the matrix
      cmatrix &a, ///< real symmetric matrix, output is replaced by orthogonal transposed matrix \e Q
      cvector1 &d, ///<diagonal elements 
      cvector2 &e ///< off-diagonal elements \e e[n-1]=0 and is irrelevant
                            )  
{
	int l,k,j,i;
//	typedef typeof(a[0][0]) cnumber;
#if __cplusplus <= 199711L
    typedef typeof(a[0][0]) cnumber;
#else
    auto xxxc=a[0][0];
    typedef decltype(xxxc) cnumber;
#endif
	cnumber scale,hh,h,g,f;
	//	int n=d.size();
    #ifdef CSTATUS_HR
    av::status S_D(n,"T");
    #ifdef FSTATUS
    S_D.SetFile();
    #endif
    #endif
	for (i=n-1;i>0;i--) {
    #ifdef CSTATUS_HR
    S_D(n-i);
    #endif
		l=i-1;
		h=scale=0.0;
		if (l > 0) {
			for (k=0;k<i;k++) scale += abs(a[k][i]);
			if (scale == 0.0)
			  e[i]=a[l][i];
			else {
				for (k=0;k<i;k++) {
				  a[k][i] /= scale;
				  	h += a[k][i]*a[k][i];
				}
				f=a[l][i];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				a[l][i]=f-g;
				f=0.0;
              //   cerr<<i<<" "<<l<<" "<<endl;
				for (j=0;j<i;j++) {
				// Next statement can be omitted if eigenvectors not wanted
				  a[i][j]=a[j][i]/h;
                 //cerr<<j<<endl;
					g=0.0;
					for (k=0;k<j+1;k++)
					  //	g += a[j][k]*a[i][k];
						g += a[k][j]*a[k][i];
					for (k=j+1;k<l+1;k++) g += a[j][k]*a[k][i];
					e[j]=g/h;
					f += e[j]*a[j][i];
					
				}
				hh=f/(h+h);
				for (j=0;j<i;j++) {
				  //f=a[i][j];
				  f=a[j][i];
					e[j]=g=e[j]-hh*f;
					for (k=0;k<j+1;k++)
					  //a[j][k] -= (f*e[k]+g*a[i][k]);
					  a[k][j] -= (f*e[k]+g*a[k][i]);
				}
			}
		} else
		  //	e[i]=a[i][l];
		e[i]=a[l][i];
		d[i]=h;
	}
	// Next statement can be omitted if eigenvectors not wanted
	d[0]=0.0;
	e[0]=0.0;
	// Contents of this loop can be omitted if eigenvectors not
	//	wanted except for statement d[i]=a[i][i];
	for (i=0;i<n;i++) {
		l=i;
		if (d[i] != 0.0) {
			for (j=0;j<l;j++) {
				g=0.0;
				for (k=0;k<l;k++)
				  //	g += a[i][k]*a[k][j];
				  	g += a[k][i]*a[j][k];
				for (k=0;k<l;k++)
				  //	a[k][j] -= g*a[k][i];
					a[j][k] -= g*a[i][k];
			}
		}
		d[i]=a[i][i];
		a[i][i]=1.0;
		for (j=0;j<l;j++) a[j][i]=a[i][j]=0.0;
	}
	for (i=1;i<n;i++) e[i-1]=e[i];  //create return so offdiagonal 0,1,..
	e[n-1]=0.0;

}


} //end tav

#endif //househoderreduction

