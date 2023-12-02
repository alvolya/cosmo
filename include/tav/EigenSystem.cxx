/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/** 
\author Alexander Volya  <http://www.volya.net>
\file EigenSystem.cxx
\author A. Volya
\note Algorithm is based on Numerical Recipes in C++, Press, et.al.
\brief Eigenvalues and Eivenvectors
\ingroup gp_nr
*/
#ifndef _EIGENSYSTEM_CXX_
#define _EIGENSYSTEM_CXX_
#include <cmath>
using std::abs;
#include "TriDiag.cxx"
#include "HouseholderReduction.cxx"
namespace tav {
/*! \brief Determine eigenvectors and eigenvalues of a real square matrix. Matrix is represented in the form
\f[ H=U D U^{T} \f]. Eigenvectors are columns of the transformation matrix
*/
 template <class cmatrix, class cvector>
 void EigenSystem(
 const int n, ///< [in] dimension of the matrix
 cmatrix &a, ///< [in,out] matrix itsef, on ouput is replaced by matrix U
 cvector &d ///< [out] vector of eigenvalues
 )
 {
// typedef typeof(a[0][0]) cnumber;
#if __cplusplus <= 199711L
     typedef typeof(a[0][0]) cnumber;
#else
     auto xxxc=a[0][0];
     typedef decltype(xxxc) cnumber;
#endif
   cnumber *e=new cnumber [n];
  tav::HouseholderReduction(n,a,d,e);
  tav::Tqli(n,d,e,a);
   delete [] e;
 }
  /*! \brief Determine eigenvectors and eigenvalues of a real square matrix. Matrix is represented in the form
\f[ H=U^{T} D U \f]. Eigenvectors are rows of the transformation matrix
*/
 template <class cmatrix, class cvector>
 void EigenSystemT(
  const int n, ///< [in] dimension of the matrix
 cmatrix &a, ///< [in,out] matrix itsef, on ouput is replaced by matrix U
 cvector &d ///< [out] vector of eigenvalues
 )
 {
//   typedef typeof(a[0][0]) cnumber;
#if __cplusplus <= 199711L
     typedef typeof(a[0][0]) cnumber;
#else
     auto xxxc=a[0][0];
     typedef decltype(xxxc) cnumber;
#endif
   cnumber *e=new cnumber [n];
  tav::HouseholderReductionT(n,a,d,e);
  tav::TqliT(n,d,e,a);
   delete [] e;
 }
} //end tav
#endif

