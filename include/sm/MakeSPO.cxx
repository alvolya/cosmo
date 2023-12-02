/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


#ifndef __MAKESPO__CXX__
#define __MAKESPO__CXX__
/*!
\author Alexander Volya  <http://www.volya.net>
\file MakeSPO.cxx
\brief Two-body interaction from system data
\todo These operators maybe combined with spherical2sp.cxx
 */
/*
This file creates a single particle matrix for two-body transitions
definition 1/4 \sum_{all indexes} V_{12;34} this V is fully antisymmetrized
10/21/2005 Write into tensor both ESP and VSP i.e. MakeESP(tensor....) 
*/
/*
Changes
10/07/2007 we use constant array aa[4] and pp[2] to assure copy in STensor
*/
 
#include <tav/tensor.cxx>
#include <tav/stensor.cxx>
using tav::Max;
using tav::Min;
#include <mav/ThreeJSymbol.cxx>
using mav::ThreeJSymbolJ;
//#include <sm/MBSQ.cxx> //depends on NextFermi...
#include <sm/sphericalint.cxx>
#include <tav/pconfiguration.cxx>
//there is a dependence on sphericalint.cxx via funciton INTS2T(V, xq, INT);
namespace SM{

///Create s.p. form of the two-body Hamiltonian
/*!
\date November 19, 2006


Operation performed is 
\f[ V_{12;34}=\sum_{L\Lambda,T\tau}\,
\sqrt{(1+\delta_{{\bf 1}{\bf 2}})(1+\delta_{{\bf 3}{\bf 4}})}\, 
V_{L\, T}^{({\bf 1}{\bf 2};{\bf 3}{\bf 4})}\, 
C_{j_{1}m_{1},j_{2}m_{2}}^{L\Lambda}C_{j_{1}m_{1},j_{2}m_{2}}^{T\tau}\, 
C_{j_{3}m_{3},j_{4}m_{4}}^{L\Lambda}C_{j_{3}m_{3},j_{4}m_{4}}^{T\tau}         
\f]   

The Hamiltonian is of the form
\f[
    H=\frac{1}{4}\sum_{0123}V_{01;23} a^\dagger_0 a^\dagger_1 a_3 a_2=
    \sum_{(01)(23)}V_{01;23} a^\dagger_0 a^\dagger_1 a_3 a_2
  \f]
- We use second summation is over pairs (01) and (23) meaning that 0<1 and 2<3.
- For T-inversal pairs are (0,1)<=(23) means that
 0<=2; if 0=2 then 1<=3; 
- Result is a map of 
spsint[4] V->first = {V[0],V[1],V[2],V[3]} -> double V->second  
- Only non-zero elements included.  

\warning The clebsh is used but phase \f$(-1)^{-2\Lambda-2\tau}\f$ 
is omitted assuming these are  integers

\todo 
- access ThreeJSymolM
- There is extra factor of 2 which comes from factor of 2 in earlier definition (ReadINT) to remove in the future!
-  Discuss 3j and Clebsh-Gordan notations

*/

int MakeVSPO(
tav::STensor<spsint,double> &VSP,  ///< Output map of elements V[0123]  
SingleParticleStates & sps, ///< s.p. states matrix
tav::Tensor<6,double> &V ///< Spherical, system data, x.V is used. 
) 
{
    //tav::STensor<spsint,double> SPV(4);
    //   VSP=0.0;
int itmp,itmpt,s,permut;
 double htmp;

 DBG(cerr<<"Generating TWO-body sp hamiltonian     "<<endl;);
//  cout<<"Generating Hamiltonian Matrix     ";
spin  Lmax=0;
spin  Tmax=0;
    for (int ll=0;ll<sps.size();ll++) {
        if (Lmax<(sps[ll].J)) Lmax=sps[ll].J;
	 if (Tmax<(sps[ll].T)) Tmax=sps[ll].T;
    }
  double LLmin,LLmax;
  double *CGCL1=new double [Lmax+1];
  double *CGCT1=new double [Tmax+1];
  double *CGCL2=new double [Lmax+1];
  double *CGCT2=new double [Tmax+1];
  int errflag;
   double TTmin,TTmax;
   int TMIN, LMIN, TMAX, LMAX, TMIN1, TMIN2, LMIN1, LMIN2; 
   //We now need to make a single particle matrix

 
   spsint aa[4];
   spsint *bb=aa+2; //just shifted reference to a[2] and a[3]
   aa[0]=0; aa[1]=1; //initialize annihilation operators
   do {
     aa[2]=aa[0];//initialize  creation 
     aa[3]=aa[1];
     do{
  
  if ( (sps[aa[0]].Jz+sps[aa[1]].Jz-sps[aa[2]].Jz-sps[aa[3]].Jz)!=0) continue;
 if ( (sps[aa[0]].Tz+sps[aa[1]].Tz-sps[aa[2]].Tz-sps[aa[3]].Tz)!=0) continue;

 
   htmp=0.0; 
  //VSP[aa[3]][aa[2]][aa[1]][aa[0]] 
	/*determine from where to where transfer (aa[2]<aa[3])<-(aa[0]<aa[1])*/
	  //Let us calculate all CGC and boundaries--------------L
	  ThreeJSymbolJ (sps[aa[2]].J/2., sps[aa[3]].J/2., sps[aa[2]].Jz/2.,sps[aa[3]].Jz/2., 
                     LLmin, LLmax, CGCL1, Lmax+1, 
	       errflag);
if (errflag!=0) continue;
	  //two permutations are used no phase trick
	  LMIN1=Round(LLmin);
	  LMAX=Round(LLmax);
	  //second one 
	  ThreeJSymbolJ (sps[aa[0]].J/2., sps[aa[1]].J/2., sps[aa[0]].Jz/2.,sps[aa[1]].Jz/2., 
                     LLmin, LLmax, CGCL2, Lmax+1, 
	       errflag);
	  if (errflag!=0) continue;
	  LMIN2=Round(LLmin);
	  LMAX=Min(Round(LLmax),LMAX); //maximum for sum
	  LMIN=Max(LMIN1,LMIN2);       //minimum for sum
	  // -------------------------------------------T
	  ThreeJSymbolJ (sps[aa[2]].T/2., sps[aa[3]].T/2., sps[aa[2]].Tz/2.,sps[aa[3]].Tz/2., 
                     TTmin, TTmax, CGCT1, Tmax+1, 
	       errflag);
if (errflag!=0) continue;


	  TMIN1=Round(TTmin);
	  TMAX=Round(TTmax);
	  //second one 
	  ThreeJSymbolJ (sps[aa[0]].T/2., sps[aa[1]].T/2., sps[aa[0]].Tz/2.,sps[aa[1]].Tz/2., 
                     TTmin, TTmax, CGCT2, Tmax+1, 
	       errflag);
if (errflag!=0) continue;
	  TMIN2=Round(TTmin);
	  TMAX=Min(Round(TTmax),TMAX); //maximum for sum
	  TMIN=Max(TMIN1,TMIN2);       //minimum for sum
	  if (LMIN>LMAX) continue;
	  if (TMIN>TMAX) continue;
	  bool addphase=
	    ((((sps[aa[1]].J-sps[aa[0]].J+sps[aa[3]].J-sps[aa[2]].J)+
	       (sps[aa[1]].T-sps[aa[0]].T+sps[aa[3]].T-sps[aa[2]].T))>>1)&1);
	  for (int l=LMIN;l<=LMAX;l++) 
	  for (int tt=TMIN;tt<=TMAX;tt++)
	    htmp+=V[sps[aa[2]].level][sps[aa[3]].level][sps[aa[0]].level][sps[aa[1]].level][l][tt]*
 CGCL1[l-LMIN1]*CGCL2[l-LMIN2]*CGCT1[tt-TMIN1]*CGCT2[tt-TMIN2]*
	   (2.0*tt+1.0)*(2.0*l+1.0)*2.0;
	
	  /*
  VSP[aa[3]][aa[2]][aa[1]][aa[0]]=htmp;   
	  VSP[aa[2]][aa[3]][aa[1]][aa[0]]=sg*htmp;  
          VSP[aa[3]][aa[2]][aa[0]][aa[1]]=sg*htmp;  
          VSP[aa[2]][aa[3]][aa[0]][aa[1]]=htmp; 
          VSP[aa[1]][aa[0]][aa[3]][aa[2]]=htmp;   
	  VSP[aa[0]][aa[1]][aa[3]][aa[2]]=sg*htmp;  
          VSP[aa[1]][aa[0]][aa[2]][aa[3]]=sg*htmp;  
          VSP[aa[0]][aa[1]][aa[2]][aa[3]]=htmp; 
	  */

	  //cout<<int(aa[0])<<" "<<int(aa[1])<<" : ";
	  //cout<<int(aa[2])<<" "<<int(aa[3])<<"  "<<htmp<<endl;
	  VSP[aa]=(addphase? -htmp:htmp);
     }while (tav::NextFermiDistribution(sps.size(), 2, (bb)));
    } while (tav::NextFermiDistribution(sps.size(), 2, aa));
 //delete [] aa;
 // cout<<VSP<<endl;

delete [] CGCL1; 
delete [] CGCT1; 
delete [] CGCL2; 
delete [] CGCT2; 

 return 1;
  };

/*! \brief Create s.p. one-body hamiltonian operator 
Assuming that \f$ H^{(I)}=\sum_{1,2} \epsilon_{1,2} a^\dagger_2 a_1 \f$ we store in 
\ref sec_cioperators matrix \f$ \epsilon_{12} \f$ which is limited to upper right or 
forward acting, namely <var>1<2</var>. This is done with the asumption of hermicity. 
\note In present spherical shell model this is irrelevant because matrix is diagonal
*/
 int MakeESPO(tav::STensor<spsint,double> &ESP,  ///< Output map of elements E[01]  
  SingleParticleStates & sps, ///< s.p. states matrix
  tav::STensor<int,double> &x ///< Spherical, system data, x.e is used. 
	     )
 {
for (spsint i=0;i<sps.size();i++) ESP(i,i)+=x(sps[i].level,sps[i].level);
return 1;	     
 } 		
  int MakeESPO(tav::STensor<spsint,double> &ESP,  ///< Output map of elements E[01]  
  SingleParticleStates & sps, ///< s.p. states matrix
  std::vector<double> &x ///< Spherical, system data, x.e is used. 
	     )
 {
for (spsint i=0;i<sps.size();i++) ESP(i,i)+=x[sps[i].level];
return 1;	     
 } 	 
 
 
 /*! \fn INT2VSPO
 \brief change tensor form of interaction to single-particle operator
*/
 int INTS2VSPO(
tav::STensor<spsint,double> &VSPX,  ///< Output map of elements V[0123]  
 tav::STensor<int,double> &INT, //two-body interaction in int form 
 SingleParticleStates & sps, ///< s.p. states matrix
 SM::ModelSpace &xq ///< model space (not neded in future)
) {
tav::Tensor<6,double> V;
INTS2T(V, xq, INT);
SM::MakeVSPO(VSPX,sps,V);
return 1;
}
  /* 

int MakeESP(
	    char *filename,   
	    SingleParticleStates & sps,
           System &x) {
   tav::Tensor<2,double> ESP(sps.size(),sps.size());
   MakeESP(ESP, sps, x);
   ESP.Write(filename);
    return 1;
}



  */
}

#endif //__MAKESPO__CXX__

