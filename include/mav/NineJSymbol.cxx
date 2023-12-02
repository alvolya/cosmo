/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/** 
\author Alexander Volya  <http://www.volya.net>
\file NineJSymbol.cxx
\brief Computation of nine-j symbols.
*/
#ifndef __NINEJSYMBOL_CXX__
#define __NINEJSYMBOL_CXX__
#include <mav/SixJSymbol.cxx>
using mav::SixJSymbol;
#include <tav/basic.cxx>
using tav::Min;
using tav::Max;
using tav::Abs;
namespace mav {
/**
9J-symbol eq. 15.19 page 133 de-Shalit, Talmi
(Messiah 1962, p. 1067; Shore and Menzel 1968, p. 282). 
  checked with a web calculator, seems fine. 
 */
double NineJSymbol(
		   double j11, double j12, double j13,
		   double j21, double j22, double j23,
		   double j31, double j32, double j33) {
  // start addressing J limits 
  int Jmin=0;
  int Jmax=10;
  Jmax=int(2.*Min(j11+j33, j21+j32, j12+j23)+0.1);
  Jmin=int(2.*Max(Abs(j11-j33), Abs(j21-j32), Abs(j12-j23))+0.1);
  //j11 j33; j21 j132; j12 j23 triangular conditions for jp
  // Jmin=2*int(fabs(j13-j24)+0.01);
  //Jmax=2*int(fabs(j13+j24)+0.01);
  //cout<<Jmin<<" mx "<<Jmax<<endl;
  double nj=0.0;
  
  //double htmp;
  double jp;
  for (int j=Jmin;j<=Jmax;j+=1) {
    jp=double(j)/2.0;
    nj+=(j+1.0)
*SixJSymbol(j11, j21, j31,
	   j32, j33, jp)
*SixJSymbol(j12, j22, j32,
	    j21, jp, j23)
*SixJSymbol(j13, j23, j33,
	    jp, j11, j12);
 
    //if (j&1) nj-=htmp;
    //else nj+=htmp;
  }
  //sum has common phase
  if (Jmin&1) return -nj;
  else return nj;
  return nj;
}



} //namespace end
//-----------------------------------------------------------------------------//
#endif 

