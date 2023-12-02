/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*!

\author Alexander Volya  <http://www.volya.net>
\file spherical2sp.cxx \brief Convert Spherical Interactions to SM::MBOperator (tav::STensor) from. 
*/
/*! \page spherical Spherical Shell Model 
\section sec_pairs Pair Operators
The pair operators are defined as 
\f[ \left(P_{L\,\Lambda,\, T\,\tau}^{({\bf 1}\,{\bf 2})}\right)^{\dagger}=\frac{1}{\sqrt{1+\delta_{{\bf 1}\,{\bf 2}}}}\,\sum_{m_{1}m_{2}t_{1}t_{2}}\, C_{j_{1}m_{1},j_{2}m_{2}}^{L\Lambda}C_{j_{1}m_{1},j_{2}m_{2}}^{T\tau}\,\, a_{1}^{\dagger}a_{2}^{\dagger}\, \f] with this definition the pair operators are normalized \f$ \langle0|P_{L\,\Lambda,\, T\,\tau}^{({\bf 1}\,{\bf 2})}\left(P_{L\,\Lambda,\, T\,\tau}^{({\bf 1}\,{\bf 2})}\right)^{\dagger}|0\rangle=1 \f$ which means that the corresponding two-particle states are normalized. The pairs are indexed by a pair of indexes denoting s.p. orbitals (shown in bold \ref sec_indexes ). Indexes can be interchanged with the following rule
\f[ \left(P_{L\,\Lambda,\, T\,\tau}^{({\bf 1}\,{\bf 2})}\right)^{\dagger}=(-)^{j_{1}+j_{2}+L+T}\,\left(P_{L\,\Lambda,\, T\,\tau}^{({\bf 2}\,{\bf 1})}\right)^{\dagger}\, \f].
\section sec_sphamiltonian Two-body Hamiltonian
The definition is
\f[  H=\sum_{L\Lambda T\tau}\sum_{({\bf 1}{\bf 2}) ({\bf 3}{\bf 4})}\,\, V_{L\, T}^{({\bf 1}{\bf 2};{\bf 3}{\bf 4})}\,\left(P_{L\,\Lambda,\, T\,\tau}^{({\bf 1}\,{\bf 2})}\right)^{\dagger}\, P_{L\,\Lambda,\, T\,\tau}^{({\bf 3}\,{\bf 4})}\,, \f]
where summation goes over pairs \f$ ({\bf 1 2}) \f$ organized that \f$ {\bf 1} \le {\bf 2} \f$
This Hamiltonian is converted to the single-particle form
\f[ V=\sum_{(12)(34)}V_{12;34}a_{1}^{\dagger}a_{2}^{\dagger}a_{4}a_{3}=\frac{1}{4}\sum_{1234}V_{12;34}a_{1}^{\dagger}a_{2}^{\dagger}a_{4}a_{3}. \f]
See discussion [Lawson 1980, p 95-96 and 27]
\todo AddPPTzVsp needs formal review, works correctly compared to oslo code
*/ 
#ifndef __SPHERICAL2SP__CXX
#define __SPHERICAL2SP__CXX
#include <tav/tensor.cxx>
#include <tav/stensor.cxx>
#include <tav/ThreeJSymbolM.cxx>
using tav::ThreeJSymbol;
#include <sm/QN.cxx>
#include <sm/SingleParticleStates.cxx>
//#include <sm/MBSQ.cxx>
#include "cioperators.cxx"
#include <tav/pconfiguration.cxx>
using SM::SPSqrM;
#include "AddPPVsp.cxx"
namespace SM{

/*! \typedef SymmetryMBOperator
\brief Define Many-body operator with additional symmetry information
*/
typedef SM::SymmetryOperator<MBOperator> SymmetryMBOperator;
  
  
/*! \brief Create s.p. form of the two-body Hamiltonian from QQ interaction 

\date September 18, 2007
- We add the the s.p. Hamiltonian a term 
\f[ 
\tilde{V}^{(1234)}_{KS}\,\sum_{\kappa \sigma} \left({\cal M}_{K\,\kappa,\, S\,\sigma}^{({\bf 1}\,{\bf 2})}\right)^{\dagger}{\cal M}_{K\,\kappa,\, S\,\sigma}^{({\bf 3}\,{\bf 4})}      \f] 
identified by parameters ijk[6]= {1,2,3,4 K,S} and V
- Multipoles are defined as 
\f[ 
({\cal M}_{K\,\kappa,\, S\,\sigma}^{({\bf 1}\,{\bf 2})})^{\dagger}=\,\sum_{m_{1}m_{2}t_{1}t_{2}}\,(-)^{j_{1}-m_{1}+1/2-t_{1}}  
\left(\begin{array}{ccc}
j_{1} & K & j_{2}\\
-m_{1} & \kappa & m_{2}\end{array}\right)\left(\begin{array}{ccc}
1/2 & S & 1/2\\
-t_{1} & \sigma & t_{2}\end{array}\right)\, a_{1}^{\dagger}a_{2}
\f] 
- output is single particle Hamiltonian given by \f$\epsilon_{12}\f$ and 
\f$V_{12;34}\f$  
The Hamiltonian is of the form
\f[
    H=\sum_{12} \epsilon_{12} a^\dagger_1 a_2 + 
    \sum_{(12),(34)}V_{12;34} a^\dagger_1 a^\dagger_2 a_4 a_3
  \f]
note that unless unless specifically defined the Hamiltonian is non-Hermitian
\see MakeSPO.cxx and function MakeVSPO()
\par Tests 
-# \f[ 
N^{({\bf 1})}=\sum_{m_{1}t_{1}}a_{1}^{\dagger}a_{1}=\sqrt{\Omega_{\bf 1}}\,{\cal M}_{0\,0,\,0\,0}^{({\bf 1}\,{\bf 1})}\quad \Omega_{\bf 1}=(2j_{\bf 1}+1)(2t_{\bf 1}+1)
\f] 
-# \f[ 
J_{\kappa}^{({\bf 1})}=\sqrt{\Omega_{\bf 1} {\bf j_1}^2}\,\left({\cal M}_{1\,\kappa,\,0\,0}^{({\bf 1}\,{\bf 1})}\right)^{\dagger}\,, \quad {\bf j}^2=j(j+1)
\f] 
The latter is in cartesian, spherical Harmonics form 
\f[ 
J^+=\sqrt{2}J_1 =\sqrt{2\Omega {\bf j}^2}\,\left({\cal M}_{1\,1,\,0\,0}\right)^{\dagger} 
\f]
In total \f[ 
\left({\cal M}_{1\,\kappa,\,0\,0}\right)^{\dagger} \, {\cal M}_{1\,\kappa,\,0\,0} =  \frac{{\bf J}^2}{\Omega {\bf j}^2} 
\f]
-# On a single-j spherical symmetry results \f$\epsilon_{12}=\delta_{12} 
\tilde{V}_K/\Omega_{\bf 1} \f$ for any multipole K.

\par Algorithm 
- Run via all s.p. states 1,2,3,4 and fill the corresponding terms based on commutation rules 
\f[ 
a^\dagger_1 a_2 a^\dagger_4 a_3 = a^\dagger_1 a_3 \delta_{24} + 
p^\dagger_{14} p_{23}\quad p_{23}=a_3 a_2.
\f] 
- order output in pairs
\todo Improve code by using ThreeJSymbolM() function
\test AddQQVsp_test.cpp
*/
int AddQQVsp(
  tav::STensor<spsint,double> &ESP,  ///< Output map of e[12] elements
  tav::STensor<spsint,double> &VSP,  ///< Output map of elements V[1234]  (12) (34)-ordered pairs
  SingleParticleStates & sps, ///< s.p. states matrix
  const int *ijk, ///< list 1 2 3 4 L T
  const double V  ///< interaction strength
              ) 
  {

    double htmp;
    bool phase;
    spsint a[4]; 
//order of labels 0123 but operators 
//a0+ a1 a3+ a2
//here we implement a brute force scan
    for (a[0]=0;a[0]<sps.size();a[0]++) {
      if (sps[a[0]].level!=ijk[0]) continue;
      for (a[1]=0;a[1]<sps.size();a[1]++) {
        if (sps[a[1]].level!=ijk[1]) continue;
        for (a[2]=0;a[2]<sps.size();a[2]++) {
          if (sps[a[2]].level!=ijk[2]) continue;
          for (a[3]=0;a[3]<sps.size();a[3]++) {
            if (sps[a[3]].level!=ijk[3]) continue;
//spin 
            if (sps[a[0]].Jz+sps[a[3]].Jz!=sps[a[1]].Jz+sps[a[2]].Jz) continue;
            if (sps[a[0]].Tz+sps[a[3]].Tz!=sps[a[1]].Tz+sps[a[2]].Tz) continue;
        
            htmp=	tav::ThreeJSymbol(sps[a[0]].J, -sps[a[0]].Jz, 
                                          2*(ijk[4]), (sps[a[0]].Jz-sps[a[1]].Jz),
                                          sps[a[1]].J, sps[a[1]].Jz);
// cout<<sps[a[0]].J/2.0 <<" "<< -sps[a[0]].Jz/2.0<<" "<< 
//     double(ijk[4])<<" "<<(sps[a[0]].Jz-sps[a[1]].Jz)/2.0<<" "<<
//     sps[a[1]].J/2.0<<" "<<sps[a[1]].Jz/2.0<<endl;
//cout<<int(a[0])<<" htmpxx="<<htmp<<" : "<<int(a[0])<<" "<<int(a[1])<<" "<<int(a[2])<<" "<<int(a[3])<<endl;
            htmp*=	tav::ThreeJSymbol(sps[a[2]].J, -sps[a[2]].Jz, 
                                          2*(ijk[4]), (sps[a[2]].Jz-sps[a[3]].Jz),
                                          sps[a[3]].J, sps[a[3]].Jz);
//cout<<int(a[0])<<" htmpx="<<htmp<<" : "<<int(a[0])<<" "<<int(a[1])<<" "<<int(a[2])<<" "<<int(a[3])<<endl;
//isospin
            htmp*=	tav::ThreeJSymbol(sps[a[0]].T, -sps[a[0]].Tz, 
                                          2*(ijk[5]), (sps[a[0]].Tz-sps[a[1]].Tz),
                                          sps[a[1]].T, sps[a[1]].Tz);
            htmp*=	tav::ThreeJSymbol(sps[a[2]].T, -sps[a[2]].Tz, 
                                          2*(ijk[5]), (sps[a[2]].Tz-sps[a[3]].Tz),
                                          sps[a[3]].T, sps[a[3]].Tz);

            if (tav::Abs(htmp)<1E-9) continue;

// cout<<int(a[0])<<" htmp="<<htmp<<" : "<<int(a[0])<<" "<<int(a[1])<<" "<<int(a[2])<<" "<<int(a[3])<<endl;
//	Pause();
//phase 

            phase=
                (((sps[a[0]].J-sps[a[0]].Jz + sps[a[2]].J-sps[a[2]].Jz + sps[a[0]].T-sps[a[0]].Tz + sps[a[2]].T-sps[a[2]].Tz)>>1)&1);   


 
            htmp=htmp*V;	    
//4 must commute via a0+ a1 a3+ a2 = a0+ a2\delta_{13} + a0+ a3+ a2 a1
            if (a[3]==a[1]) ESP(a[0],a[2])+=(phase ? -htmp : htmp);	

//Fermi forbidden pairs
            if (a[0]==a[3]) continue;
            if (a[1]==a[2]) continue; 
            spsint b[4];
            b[0]=a[0];b[1]=a[3]; b[2]=a[1]; b[3]=a[2]; 
 //note extra switch in operator a0+ a3+ a2 a1 = p{02}+ p{12} 
 //the pair-state is b[0123]<-a[0312]
 // the order is now a0+a3+ a2 a1, note that last two are switched so no minus
//create pairs by ordering operators
            if(b[0]>b[1]) {tav::Swap(b[0],b[1]); phase=!phase;}
            if(b[2]>b[3]) {tav::Swap(b[2],b[3]); phase=!phase;}
// Uncomment below to order pairs (01)<(23), assuming hermicity
/*
            if (b[0]>b[2]) {tav::Swap(b[0],b[2]); tav::Swap(b[1],b[3]);}
            if ((b[0]==b[2])&&(b[1]>b[3])) {tav::Swap(b[0],b[2]); tav::Swap(b[1],b[3]);}
*/
//enter data 
            VSP[b]+=(phase ? -htmp : htmp); 
   
          }}}} 
          return 1;
  } //end of function 
  
/*! \brief Compatibility form, add s.p. in Tensor form
  */
int AddQQVsp(
tav::Tensor<2,double> &ESP,  ///< Output double-array of e[1][2]
tav::STensor<spsint,double> &VSP,  ///< Output map of elements V[1234]  (12) (34)-ordered pairs
SingleParticleStates & sps, ///< s.p. states matrix
const int *ijk, ///< list 1 2 3 4 L T
const double V  ///< interaction strength
) 
{
  

    double htmp;
    bool phase;
    spsint a[4]; 
//order of labels 0123 but operators 
//a0+ a1 a3+ a2
//here we implement a brute force scan
    for (a[0]=0;a[0]<sps.size();a[0]++) {
	if (sps[a[0]].level!=ijk[0]) continue;
    for (a[1]=0;a[1]<sps.size();a[1]++) {
	if (sps[a[1]].level!=ijk[1]) continue;
    for (a[2]=0;a[2]<sps.size();a[2]++) {
	if (sps[a[2]].level!=ijk[2]) continue;
    for (a[3]=0;a[3]<sps.size();a[3]++) {
	if (sps[a[3]].level!=ijk[3]) continue;
//spin 
	if (sps[a[0]].Jz+sps[a[3]].Jz!=sps[a[1]].Jz+sps[a[2]].Jz) continue;
	if (sps[a[0]].Tz+sps[a[3]].Tz!=sps[a[1]].Tz+sps[a[2]].Tz) continue;
        
htmp=	tav::ThreeJSymbol(sps[a[0]].J, -sps[a[0]].Jz, 
		  2*(ijk[4]), (sps[a[0]].Jz-sps[a[1]].Jz),
			  sps[a[1]].J, sps[a[1]].Jz);
// cout<<sps[a[0]].J/2.0 <<" "<< -sps[a[0]].Jz/2.0<<" "<< 
//     double(ijk[4])<<" "<<(sps[a[0]].Jz-sps[a[1]].Jz)/2.0<<" "<<
//     sps[a[1]].J/2.0<<" "<<sps[a[1]].Jz/2.0<<endl;
//cout<<int(a[0])<<" htmpxx="<<htmp<<" : "<<int(a[0])<<" "<<int(a[1])<<" "<<int(a[2])<<" "<<int(a[3])<<endl;
htmp*=	tav::ThreeJSymbol(sps[a[2]].J, -sps[a[2]].Jz, 
		  2*(ijk[4]), (sps[a[2]].Jz-sps[a[3]].Jz),
			  sps[a[3]].J, sps[a[3]].Jz);
//cout<<int(a[0])<<" htmpx="<<htmp<<" : "<<int(a[0])<<" "<<int(a[1])<<" "<<int(a[2])<<" "<<int(a[3])<<endl;
//isospin
htmp*=	tav::ThreeJSymbol(sps[a[0]].T, -sps[a[0]].Tz, 
		  2*(ijk[5]), (sps[a[0]].Tz-sps[a[1]].Tz),
			  sps[a[1]].T, sps[a[1]].Tz);
htmp*=	tav::ThreeJSymbol(sps[a[2]].T, -sps[a[2]].Tz, 
		  2*(ijk[5]), (sps[a[2]].Tz-sps[a[3]].Tz),
			  sps[a[3]].T, sps[a[3]].Tz);

 if (tav::Abs(htmp)<1E-9) continue;

// cout<<int(a[0])<<" htmp="<<htmp<<" : "<<int(a[0])<<" "<<int(a[1])<<" "<<int(a[2])<<" "<<int(a[3])<<endl;
//	Pause();
//phase 

 phase=
     (((sps[a[0]].J-sps[a[0]].Jz + sps[a[2]].J-sps[a[2]].Jz + sps[a[0]].T-sps[a[0]].Tz + sps[a[2]].T-sps[a[2]].Tz)>>1)&1);   


 
 htmp=htmp*V;	    
//4 must commute via a0+ a1 a3+ a2 = a0+ a2\delta_{13} + a0+ a3+ a2 a1
 if (a[3]==a[1]) ESP[a[0]][a[2]]+=(phase ? -htmp : htmp);	

//Fermi forbidden pairs
 if (a[0]==a[3]) continue;
 if (a[1]==a[2]) continue; 
 spsint b[4];
 b[0]=a[0];b[1]=a[3]; b[2]=a[1]; b[3]=a[2]; 
 //note extra switch in operator a0+ a3+ a2 a1 = p{02}+ p{12} 
 //the pair-state is b[0123]<-a[0312]
 // the order is now a0+a3+ a2 a1, note that last two are switched so no minus
//create pairs by ordering operators
 if(b[0]>b[1]) {tav::Swap(b[0],b[1]); phase=!phase;}
 if(b[2]>b[3]) {tav::Swap(b[2],b[3]); phase=!phase;}
// Uncomment below to order pairs (01)<(23), assuming hermicity
/*
 if (b[0]>b[2]) {tav::Swap(b[0],b[2]); tav::Swap(b[1],b[3]);}
 if ((b[0]==b[2])&&(b[1]>b[3])) {tav::Swap(b[0],b[2]); tav::Swap(b[1],b[3]);}
*/
//enter data 
 VSP[b]+=(phase ? -htmp : htmp); 
   
    }}}} 
    return 1;
} //end of function 

/*! \brief Create a single-particle operator from a spherical interaction, \f$ V_{34;12} \f$ 
The part of two-body interaction that depends on spherical levels \f$ {\bf 12;34} \f$.
We enforce forward ordering here \f$ {\bf (12)\le(34) } \f$
\f[ \acute{H}_{LT}({\bf 34;12})=\frac{V_{LT}({\bf 34;12})}{4} (1+\delta_{\bf 1 2})(1+\delta_{\bf 3 4}) \sum_{M \tau} P^\dagger_{LM T\tau}({\bf 3,4}) P_{LM T\tau}({\bf 1,2})= \frac{1}{4}\sum_{1234} V_{34;12}  a^\dagger_3 a^\dagger_4 a_2 a_1 
\f]
\f[ V_{34;12}=V_{LT}({\bf 12;34}) \sqrt{(1+\delta_{\bf 1 2})(1+\delta_{\bf 3 4})} C^{LT}_{3,4} C^{LT}_{12}  \f] 
encoded as <var> {1234} </var>. Operator forward ordered.
\f[ \left(P_{L\,\Lambda,\, T\,\tau}^{({\bf 1}\,{\bf 2})}\right)^{\dagger}=(-)^{j_{1}+j_{2}+L+T}\,\left(P_{L\,\Lambda,\, T\,\tau}^{({\bf 2}\,{\bf 1})}\right)^{\dagger} \f]
\todo This is not the most effecient and elegant code, Use full list of ClebshGordan
pair ordering..change forward convention? 

\example SphericalPPINT_test.cpp
boson capable
This is original version, more transparent and robust
 */
int AddPPVsp1(
tav::STensor<spsint,double> &VSP,  ///< Output map of elements V[1234]  (12) (34)-ordered pairs
SingleParticleStates & sps, ///< s.p. states matrix
const int *ijk, ///< list 1 2 3 4 L T
const double V  ///< interaction strength
) 
{
    double htmp;
    bool phase;
    //check ordering of input operator
    spsint a[4];  spsint b[4]; 
//order of labels 0123 but operators 
//a3+ a2+ a1 a0
//here we implement a brute force scan
    for (a[0]=0;a[0]<sps.size();a[0]++) {
	if (sps[a[0]].level!=ijk[0]) continue;
    for (a[1]=0;a[1]<sps.size();a[1]++) {
	if (sps[a[1]].level!=ijk[1]) continue;
        if (((sps[a[0]].J)&1)&&(a[1]==a[0])) continue; //boson fermion
    for (a[2]=0;a[2]<sps.size();a[2]++) { 
//pairs are NOT forward ordered: we can not do that because of mapping to levels
	if (sps[a[2]].level!=ijk[2]) continue;
    for (a[3]=0;a[3]<sps.size();a[3]++) {
	if (sps[a[3]].level!=ijk[3]) continue;
        if (((sps[a[2]].J)&1)&&(a[2]==a[3])) continue; //fermion and boson
//spin        
	if (sps[a[0]].Jz+sps[a[1]].Jz!=sps[a[3]].Jz+sps[a[2]].Jz) continue;
	if (sps[a[0]].Tz+sps[a[1]].Tz!=sps[a[3]].Tz+sps[a[2]].Tz) continue;
        //cout<<int(a[0])<<" "<<int(a[1])<<" "<<int(a[2])<<" "<<int(a[3])<<endl;
        phase=false;
        htmp=tav::ClebschGordan(sps[a[0]].J, sps[a[0]].Jz, 
                                   sps[a[1]].J, sps[a[1]].Jz,
		                   2*(ijk[4]), (sps[a[0]].Jz+sps[a[1]].Jz));
// cout<<sps[a[0]].J/2.0 <<" "<< -sps[a[0]].Jz/2.0<<" "<< 
//     double(ijk[4])<<" "<<(sps[a[0]].Jz-sps[a[1]].Jz)/2.0<<" "<<
//     sps[a[1]].J/2.0<<" "<<sps[a[1]].Jz/2.0<<endl;
//cout<<int(a[0])<<" htmpxx="<<htmp<<" : "<<int(a[0])<<" "<<int(a[1])<<" "<<int(a[2])<<" "<<int(a[3])<<endl;
htmp*=	tav::ClebschGordan(sps[a[2]].J, sps[a[2]].Jz, 
                                   sps[a[3]].J, sps[a[3]].Jz,
		                   2*(ijk[4]), (sps[a[2]].Jz+sps[a[3]].Jz));
//cout<<int(a[0])<<" htmpx="<<htmp<<" : "<<int(a[0])<<" "<<int(a[1])<<" "<<int(a[2])<<" "<<int(a[3])<<endl;
//isospin
 htmp*=tav::ClebschGordan(sps[a[0]].T, sps[a[0]].Tz, 
                         sps[a[1]].T, sps[a[1]].Tz,
                             2*(ijk[5]), (sps[a[0]].Tz+sps[a[1]].Tz));
 htmp*=	tav::ClebschGordan(sps[a[2]].T, sps[a[2]].Tz, 
                           sps[a[3]].T, sps[a[3]].Tz,
                               2*(ijk[5]), (sps[a[2]].Tz+sps[a[3]].Tz));

 if (tav::Abs(htmp)<1E-9) continue;

 //cout<<htmp<<endl;
 htmp=htmp*V/sqrt((double((ijk[0]==ijk[1]))+1.)*(double((ijk[2]==ijk[3]))+1.));
 /* 
 Here we use interaction in isospin format, V is the same but levels are indexed separately
 For interaction in the pn format if I have a pn pair from the same level but they are indexed separately then 
 we need to multiply by sqrt(2)
 */
 if ((ijk[0]!=ijk[1])&&(sps[a[0]].N==sps[a[1]].N)&&(sps[a[0]].J==sps[a[1]].J)&&(sps[a[0]].Tz==-sps[a[1]].Tz)&&(sps[a[0]].P==sps[a[1]].P)) htmp*=sqrt(2.);
 if ((ijk[2]!=ijk[3])&&(sps[a[2]].N==sps[a[3]].N)&&(sps[a[2]].J==sps[a[3]].J)&&(sps[a[2]].Tz==-sps[a[3]].Tz)&&(sps[a[2]].P==sps[a[3]].P)) htmp*=sqrt(2.);
 //end scaling
b[0]=a[0];b[1]=a[1];b[2]=a[2];b[3]=a[3]; 
// if(b[0]>b[1]) {tav::Swap(b[0],b[1]); phase=!phase;}
// if(b[2]>b[3]) {tav::Swap(b[2],b[3]); phase=!phase;}
if(b[0]>b[1]) {tav::Swap(b[0],b[1]); phase  ^=(sps[a[0]].J)&1;}
if(b[2]>b[3]) {tav::Swap(b[2],b[3]);  phase ^=(sps[a[2]].J)&1;} //check for boson phase
 
  //set forward ordering of pairs; use Array comparison for that
 tav::ArrayComparison<spsint> X(2);
 if (X(b+2,b)) {tav::Swap(b[0],b[2]); tav::Swap(b[1],b[3]);}
 if ((ijk[0]==ijk[2])&&(ijk[1]==ijk[3])) htmp/=2.0;
 if ((b[0]==b[2])&&(b[1]==b[3])) htmp*=2.0; //my bad forward convention
 ///\todo change forward convention?
 VSP[b]+=(phase ? -htmp : htmp); 
 
    }}}} 
    return 1;
} //end of function

/*! \brief  add interaction matrix but require pp, pn or nn only
boson capable
*/
int AddPPVsp(
tav::STensor<spsint,double> &VSP,  ///< Output map of elements V[1234]  (12) (34)-ordered pairs
SingleParticleStates & sps, ///< s.p. states matrix
const int *ijk, ///< list 1 2 3 4 L T
const double V,  ///< interaction strength
int pnstate ///<extra condition on pn-state 1-NN, 0-PN, -1 PP
) 
{
    double htmp;
    bool phase;
    //check ordering of input operator
    spsint a[4];  spsint b[4]; 
//order of labels 0123 but operators 
//a3+ a2+ a1 a0
//here we implement a brute force scan
    for (a[0]=0;a[0]<sps.size();a[0]++) {
	if (sps[a[0]].level!=ijk[0]) continue;
    for (a[1]=0;a[1]<sps.size();a[1]++) {
	if (sps[a[1]].level!=ijk[1]) continue;
        //if (a[1]==a[0]) continue;
        if (((sps[a[0]].J)&1)&&(a[1]==a[0])) continue; //boson fermion
    for (a[2]=0;a[2]<sps.size();a[2]++) { 
//pairs are NOT forward ordered: we can not do that because of mapping to levels
	if (sps[a[2]].level!=ijk[2]) continue;
    for (a[3]=0;a[3]<sps.size();a[3]++) {
	if (sps[a[3]].level!=ijk[3]) continue;
        //if (a[2]==a[3]) continue;
        if (((sps[a[2]].J)&1)&&(a[2]==a[3])) continue; //fermion and boson
//spin        
	if (sps[a[0]].Jz+sps[a[1]].Jz!=sps[a[3]].Jz+sps[a[2]].Jz) continue;
	if (sps[a[0]].Tz+sps[a[1]].Tz!=sps[a[3]].Tz+sps[a[2]].Tz) continue;
    if (int(sps[a[0]].Tz+sps[a[1]].Tz)!=2*pnstate) continue; //this is extra condition on pn state
        //cout<<int(a[0])<<" "<<int(a[1])<<" "<<int(a[2])<<" "<<int(a[3])<<endl;
        phase=false;
        htmp=tav::ClebschGordan(sps[a[0]].J, sps[a[0]].Jz, 
                                   sps[a[1]].J, sps[a[1]].Jz,
		                   2*(ijk[4]), (sps[a[0]].Jz+sps[a[1]].Jz));
// cout<<sps[a[0]].J/2.0 <<" "<< -sps[a[0]].Jz/2.0<<" "<< 
//     double(ijk[4])<<" "<<(sps[a[0]].Jz-sps[a[1]].Jz)/2.0<<" "<<
//     sps[a[1]].J/2.0<<" "<<sps[a[1]].Jz/2.0<<endl;
//cout<<int(a[0])<<" htmpxx="<<htmp<<" : "<<int(a[0])<<" "<<int(a[1])<<" "<<int(a[2])<<" "<<int(a[3])<<endl;
htmp*=	tav::ClebschGordan(sps[a[2]].J, sps[a[2]].Jz, 
                                   sps[a[3]].J, sps[a[3]].Jz,
		                   2*(ijk[4]), (sps[a[2]].Jz+sps[a[3]].Jz));
//cout<<int(a[0])<<" htmpx="<<htmp<<" : "<<int(a[0])<<" "<<int(a[1])<<" "<<int(a[2])<<" "<<int(a[3])<<endl;
//isospin
 htmp*=tav::ClebschGordan(sps[a[0]].T, sps[a[0]].Tz, 
                         sps[a[1]].T, sps[a[1]].Tz,
                             2*(ijk[5]), (sps[a[0]].Tz+sps[a[1]].Tz));
 htmp*=	tav::ClebschGordan(sps[a[2]].T, sps[a[2]].Tz, 
                           sps[a[3]].T, sps[a[3]].Tz,
                               2*(ijk[5]), (sps[a[2]].Tz+sps[a[3]].Tz));

 if (tav::Abs(htmp)<1E-9) continue;

 //cout<<htmp<<endl;
 htmp=htmp*V/sqrt((double((ijk[0]==ijk[1]))+1.)*(double((ijk[2]==ijk[3]))+1.));
 /* 
 Here we use interaction in isospin format, V is the same but levels are indexed separately
 For interaction in the pn format if I have a pn pair from the same level but they are indexed separately then 
 we need to multiply by sqrt(2)
 */
 if ((ijk[0]!=ijk[1])&&(sps[a[0]].N==sps[a[1]].N)&&(sps[a[0]].J==sps[a[1]].J)&&(sps[a[0]].Tz==-sps[a[1]].Tz)&&(sps[a[0]].P==sps[a[1]].P)) htmp*=sqrt(2.);
 if ((ijk[2]!=ijk[3])&&(sps[a[2]].N==sps[a[3]].N)&&(sps[a[2]].J==sps[a[3]].J)&&(sps[a[2]].Tz==-sps[a[3]].Tz)&&(sps[a[2]].P==sps[a[3]].P)) htmp*=sqrt(2.);
 //end scaling
b[0]=a[0];b[1]=a[1];b[2]=a[2];b[3]=a[3]; 
 //order matrix elements and swap phase based on CGC
 // if(b[0]>b[1]) {tav::Swap(b[0],b[1]); phase=!phase;}
// if(b[2]>b[3]) {tav::Swap(b[2],b[3]); phase=!phase;}
 if(b[0]>b[1]) {tav::Swap(b[0],b[1]); phase  ^=(sps[a[0]].J)&1;}
 if(b[2]>b[3]) {tav::Swap(b[2],b[3]);  phase ^=(sps[a[2]].J)&1;} //check for boson phase

  //set forward ordering of pairs; use Array comparison for that
 tav::ArrayComparison<spsint> X(2);
 if (X(b+2,b)) {tav::Swap(b[0],b[2]); tav::Swap(b[1],b[3]);}
 if ((ijk[0]==ijk[2])&&(ijk[1]==ijk[3])) htmp/=2.0;
 if ((b[0]==b[2])&&(b[1]==b[3])) htmp*=2.0; //my bad forward convention
 ///\todo change forward convention?
 VSP[b]+=(phase ? -htmp : htmp); 
 
    }}}} 
    return 1;
} //end of function 



/*! \brief Create a single-particle operator from a spherical interaction, \f$ V_{34;12} \f$ 
but with fixed Tz projections, so separate interaction for nn, np and pp.
 */
int AddPPTzVsp(
tav::STensor<spsint,double> &VSP,  ///< Output map of elements V[1234]  (12) (34)-ordered pairs
SingleParticleStates & sps, ///< s.p. states matrix
const int *ijk, ///< list 1 2 3 4 L Tz
const double V  ///< interaction strength
) 
{
    double htmp;
    bool phase;
    //check ordering of input operator
    spsint a[4];  spsint b[4]; 
//order of labels 0123 but operators 
//a3+ a2+ a1 a0
//here we implement a brute force scan
    for (a[0]=0;a[0]<sps.size();a[0]++) {
	if (sps[a[0]].level!=ijk[0]) continue;
    for (a[1]=0;a[1]<sps.size();a[1]++) {
	if (sps[a[1]].level!=ijk[1]) continue;
        if (a[1]==a[0]) continue;
    for (a[2]=0;a[2]<sps.size();a[2]++) { 
//pairs are NOT forward ordered: we can not do that because of mapping to levels
	if (sps[a[2]].level!=ijk[2]) continue;
    for (a[3]=0;a[3]<sps.size();a[3]++) {
	if (sps[a[3]].level!=ijk[3]) continue;
        if (a[2]==a[3]) continue;
//spin        
	if (sps[a[0]].Jz+sps[a[1]].Jz!=sps[a[3]].Jz+sps[a[2]].Jz) continue;
	if (sps[a[0]].Tz+sps[a[1]].Tz!=sps[a[3]].Tz+sps[a[2]].Tz) continue;
	if (sps[a[0]].Tz+sps[a[1]].Tz!=2*ijk[5]) continue; //note sign of isospin operator
        //cout<<int(a[0])<<" "<<int(a[1])<<" "<<int(a[2])<<" "<<int(a[3])<<endl;
        phase=false;
        htmp=tav::ClebschGordan(sps[a[0]].J, sps[a[0]].Jz, 
                                   sps[a[1]].J, sps[a[1]].Jz,
		                   2*(ijk[4]), (sps[a[0]].Jz+sps[a[1]].Jz));
// cout<<sps[a[0]].J/2.0 <<" "<< -sps[a[0]].Jz/2.0<<" "<< 
//     double(ijk[4])<<" "<<(sps[a[0]].Jz-sps[a[1]].Jz)/2.0<<" "<<
//     sps[a[1]].J/2.0<<" "<<sps[a[1]].Jz/2.0<<endl;
//cout<<int(a[0])<<" htmpxx="<<htmp<<" : "<<int(a[0])<<" "<<int(a[1])<<" "<<int(a[2])<<" "<<int(a[3])<<endl;
htmp*=	tav::ClebschGordan(sps[a[2]].J, sps[a[2]].Jz, 
                                   sps[a[3]].J, sps[a[3]].Jz,
		                   2*(ijk[4]), (sps[a[2]].Jz+sps[a[3]].Jz));
//cout<<int(a[0])<<" htmpx="<<htmp<<" : "<<int(a[0])<<" "<<int(a[1])<<" "<<int(a[2])<<" "<<int(a[3])<<endl;
//isospin

 if (tav::Abs(htmp)<1E-9) continue;

 //cout<<htmp<<endl;
 
 
 htmp*=V;
 //There is some guess work here, if pn, no pair normalization
 if (ijk[5]!=0) htmp/=sqrt((double((ijk[0]==ijk[1]))+1.)*(double((ijk[2]==ijk[3]))+1.));
 
 
 //FatalError("pn part is not working correctly identical fermion problem");
 /*to deal with isospin we isert (pn)=P(T=0)+P(T=1) which means with substituted cgc multiplication by
 (tz_1-tz_2+1)*(tz_3-tz_4+1)/2.0 only when total pair tz=0 another 1/2 from pn normalization
 */
 if (ijk[5]==0) htmp*=0.25*((sps[a[0]].Tz-sps[a[1]].Tz)/2.0+1.0)*((sps[a[2]].Tz-sps[a[3]].Tz)/2.0+1.0); 
b[0]=a[0];b[1]=a[1];b[2]=a[2];b[3]=a[3]; 
 if(b[0]>b[1]) {tav::Swap(b[0],b[1]); phase=!phase;}
 if(b[2]>b[3]) {tav::Swap(b[2],b[3]); phase=!phase;}
//phase^=((sps[a[2]].J+sps[a[3]].J)/2-ijk[4])&1; //modify phase to match cgc in angular momentum
  //set forward ordering of pairs; use Array comparison for that
 tav::ArrayComparison<spsint> X(2);
 if (X(b+2,b)) {tav::Swap(b[0],b[2]); tav::Swap(b[1],b[3]);}
 //if (ijk[5]==0) htmp/=2.0;
 //else 
 if ((ijk[0]==ijk[2])&&(ijk[1]==ijk[3])) htmp/=2.0;
 if ((b[0]==b[2])&&(b[1]==b[3])) htmp*=2.0; //my bad forward convention
 ///\todo change forward convention?
 VSP[b]+=(phase ? -htmp : htmp); 
 //cerr<<"matrix element "<<V<<" "<<htmp<<endl;
    }}}} 
    return 1;
} //end of function 



/*! \brief Add H(34,12) but with fixed projection of P_LM operator.
See SM::AddPPVsp(), everything is the same but we add a particular projection 
\f[ \acute{H}_{LM\,T}({\bf 34;12})=\frac{V_{LM\,T}({\bf 34;12})}{4} (1+\delta_{\bf 1 2})(1+\delta_{\bf 3 4}) \sum_{\tau} P^\dagger_{LM T\tau}({\bf 3,4}) P_{LM T\tau}({\bf 1,2})= \frac{1}{4}\sum_{1234} V_{34;12}  a^\dagger_3 a^\dagger_4 a_2 a_1 
\f]
\par Purpose
The hermitian operator \f$ {H}_M({\bf 34;12}) \f$ does not conserve angular momentum. However expectation value in  spin-zero states \f$ \langle \alpha|H| \beta \rangle=(2L+1) \langle \alpha|H_{M}| \beta \rangle\f$. Thus with large parameter \f$ \lambda \f$
spin zero eginvalues of \f$ (2L+1)H_M+ \lambda {\bf J}^2 \f$ are equal to eigenvalues
of the original Hamiltonian. 
*/
int AddPPVMsp(int M, ///< we add P_LM  component with given M
    tav::STensor<spsint,double> &VSP,  ///< Output map of elements V[1234]  (12) (34)-ordered pairs
SingleParticleStates & sps, ///< s.p. states matrix
const int *ijk, ///< list 1 2 3 4 L T
const double V  ///< interaction strength
            ) 
{

  double htmp;
  bool phase;
  spsint a[4]; 
//order of labels 0123 but operators 
//a3+ a2+ a1 a0
//here we implement a brute force scan
  for (a[0]=0;a[0]<sps.size();a[0]++) {
    if (sps[a[0]].level!=ijk[0]) continue;
    for (a[1]=a[0]+1;a[1]<sps.size();a[1]++) {
      if (sps[a[1]].level!=ijk[1]) continue;
     //pair (12) is ready at this stage, force projection L
      if (sps[a[0]].Jz+sps[a[1]].Jz!=2*M) continue;
      for (a[2]=a[0];a[2]<sps.size();a[2]++) { //pairs are forward ordered
        if (sps[a[2]].level!=ijk[2]) continue;
        for (a[3]=((a[0]==a[2])? a[1]: a[2]+1);a[3]<sps.size();a[3]++) {
          if (sps[a[3]].level!=ijk[3]) continue;
          if (sps[a[2]].Jz+sps[a[3]].Jz!=2*M) continue;
//spin not needed
          //if (sps[a[0]].Jz+sps[a[1]].Jz!=sps[a[3]].Jz+sps[a[2]].Jz) continue;
          if (sps[a[0]].Tz+sps[a[1]].Tz!=sps[a[3]].Tz+sps[a[2]].Tz) continue;
          phase=false;
          htmp=tav::ClebschGordan(sps[a[0]].J, sps[a[0]].Jz, 
                                  sps[a[1]].J, sps[a[1]].Jz,
                                  2*(ijk[4]), (sps[a[0]].Jz+sps[a[1]].Jz));
// cout<<sps[a[0]].J/2.0 <<" "<< -sps[a[0]].Jz/2.0<<" "<< 
//     double(ijk[4])<<" "<<(sps[a[0]].Jz-sps[a[1]].Jz)/2.0<<" "<<
//     sps[a[1]].J/2.0<<" "<<sps[a[1]].Jz/2.0<<endl;
//cout<<int(a[0])<<" htmpxx="<<htmp<<" : "<<int(a[0])<<" "<<int(a[1])<<" "<<int(a[2])<<" "<<int(a[3])<<endl;
          htmp*=	tav::ClebschGordan(sps[a[2]].J, sps[a[2]].Jz, 
                                           sps[a[3]].J, sps[a[3]].Jz,
                                           2*(ijk[4]), (sps[a[2]].Jz+sps[a[3]].Jz));
//cout<<int(a[0])<<" htmpx="<<htmp<<" : "<<int(a[0])<<" "<<int(a[1])<<" "<<int(a[2])<<" "<<int(a[3])<<endl;
//isospin
          htmp*=tav::ClebschGordan(sps[a[0]].T, sps[a[0]].Tz, 
                                   sps[a[1]].T, sps[a[1]].Tz,
                                   2*(ijk[5]), (sps[a[0]].Tz+sps[a[1]].Tz));
          htmp*=	tav::ClebschGordan(sps[a[2]].T, sps[a[2]].Tz, 
                                           sps[a[3]].T, sps[a[3]].Tz,
                                           2*(ijk[5]), (sps[a[2]].Tz+sps[a[3]].Tz));

          if (tav::Abs(htmp)<1E-9) continue;

 //cout<<htmp<<endl;
          htmp=htmp*V*sqrt((double((ijk[0]==ijk[1]))+1.)*(double((ijk[2]==ijk[3]))+1.));
 //note we do not divide by 4 because pairs are commuted to my order 
 //note extra switch in operator a0+ a3+ a2 a1 = p{02}+ p{12} 
 //the pair-state is b[0123]<-a[0312]
 // the order is now a0+a3+ a2 a1, note that last two are switched so no minus
//create pairs by ordering operators not needed in this case
// if(b[0]>b[1]) {tav::Swap(b[0],b[1]); phase=!phase;}
// if(b[2]>b[3]) {tav::Swap(b[2],b[3]); phase=!phase;}
// Uncomment below to order pairs (01)<(23), assuming hermicity
/*
          if (b[0]>b[2]) {tav::Swap(b[0],b[2]); tav::Swap(b[1],b[3]);}
          if ((b[0]==b[2])&&(b[1]>b[3])) {tav::Swap(b[0],b[2]); tav::Swap(b[1],b[3]);}
*/
//enter data 
// cerr<<int(a[0])<<" "<<int(a[1])<<" "<<int(a[2])<<" "<<int(a[3])<<" "<<htmp<<endl;
         // 
          VSP[a]+=(phase ? -htmp : htmp); 
 
        }}}} 
        return 1;
} //end of function 
/*! \brief Generate \f$J^+\f$ operator

The single particle operator 
\f$J^+=\sum_{jm} a^\dagger_{jm+1} a_{jm} \sqrt{(j-m)(j+m+1)} = \sum_{12}q_{12} a^\dagger_2 a_1\f$
The return is a set of \f$q_{12}=\sqrt{(j_1-m_1)(j_1+m_1+1)} \f$ in the form of STensor, where \f$m_2=m_1+1\f$ so that 
 
 */
  /*this program generates two-body hamiltonian for J+ operator
  remember that first goes is annihilation operator, then creation
  we have \sum_jm a^+(jm+1) a(jm) \sqrt{(j-m)(j+m+1)}
  08/07/08 We fix problem conjugated order in computation
  */
  int SPJPlus(tav::STensor<spsint,double> &JO,   
	    SingleParticleStates & sps) {
    JO.clear();
    spsint pp[2]; //we first do pairs, must use constants
    pp[0]=0; pp[1]=1; //initialize pairs 
    //m+1 is always futher than m so everything is correct despite that this is
    //a_jm and a_jm-1 generation
    do {
      if (sps[pp[0]].level!=sps[pp[1]].level) continue; //same level
      if (sps[pp[0]].Tz!=sps[pp[1]].Tz) continue; //same prijection Tz
      if ((sps[pp[0]].Jz+2)==sps[pp[1]].Jz) 
	{
	  //this is a normal order operator p[0]=jm appears first 
JO[pp]= sqrt((sps[pp[0]].J-sps[pp[0]].Jz)*(sps[pp[0]].J+sps[pp[0]].Jz+2.))/2.0;
	  continue;}
      if ((sps[pp[0]].Jz-2)==sps[pp[1]].Jz) 
	{//this wrong order pp[1]=jm 
          tav::Swap(pp[0],pp[1]);
JO[pp]=sqrt((sps[pp[0]].J-sps[pp[0]].Jz)*(sps[pp[0]].J+sps[pp[0]].Jz+2.))/2.0;
 tav::Swap(pp[0],pp[1]);
	  continue;}

    } while (tav::NextFermiDistribution(sps.size(), 2, pp));
   
      //    Pause();

 //delete [] pp;


 return 1;
  };


  /*!T+ operator, same as J+ with J<->T change
   */
  int SPTPlus(tav::STensor<spsint,double> &JO,   
	    SingleParticleStates & sps) {
    JO.clear();
    spsint pp[2]; //we first do pairs
    pp[0]=0; pp[1]=1; //initialize pairs 
    //m+1 is always futher than m so everything is correct despite that this is
    //a_jm and a_jm-1 generation
    do {
     // if (sps[pp[0]].level!=sps[pp[1]].level) continue; //same level
      if (sps[pp[0]].Jz!=sps[pp[1]].Jz) continue; //same prijection Tz
      if (sps[pp[0]].N!=sps[pp[1]].N) continue;
       if (sps[pp[0]].P!=sps[pp[1]].P) continue;
       if (sps[pp[0]].L!=sps[pp[1]].L) continue;
        if (sps[pp[0]].T!=sps[pp[1]].T) continue;
         if (sps[pp[0]].J!=sps[pp[1]].J) continue;
      if ((sps[pp[0]].Tz+2)==sps[pp[1]].Tz) 
	{
	  //this is a normal order operator p[0]=jm appears first 
JO[pp]= sqrt((sps[pp[0]].T-sps[pp[0]].Tz)*(sps[pp[0]].T+sps[pp[0]].Tz+2.))/2.0;
	  continue;}
      if ((sps[pp[0]].Tz-2)==sps[pp[1]].Tz) 
	{//this wrong order pp[1]=jm 
          tav::Swap(pp[0],pp[1]);
JO[pp]=sqrt((sps[pp[0]].T-sps[pp[0]].Tz)*(sps[pp[0]].T+sps[pp[0]].Tz+2.))/2.0;
tav::Swap(pp[0],pp[1]);
//cout<<JO[pp]<<" "<<int(pp[0])<<" "<<int(pp[1])<<endl;	  
continue;}

    } while (tav::NextFermiDistribution(sps.size(), 2, pp));
   
      //    Pause();



 return 1;
  };
		

  int SPJJ(tav::STensor<spsint,double> &JSP,  tav::STensor<spsint,double> &VSP,  
	    SingleParticleStates & sps) {

    tav::STensor<spsint,double> JO(2);
  SM::SPJPlus(JO, sps);
 SM::SPSqrM(JSP,VSP,JO);
SM::StripHermitianOperator(JSP);
 SM::StripHermitianOperator(VSP);
 return 1;
  };

  int SPTT(tav::STensor<spsint,double> &JSP,  tav::STensor<spsint,double> &VSP,  
	    SingleParticleStates & sps) {

    tav::STensor<spsint,double> JO(2);
  SM::SPTPlus(JO, sps);
 SM::SPSqrM(JSP,VSP,JO);
SM::StripHermitianOperator(JSP);
 SM::StripHermitianOperator(VSP);
 return 1;
  };


/*! \brief Add/create multipole operator.  
Definition: 
\f[ 
 ({\cal M}_{K\,\kappa,\, S\,\sigma}^{({\bf 1}\,{\bf 2})})^{\dagger}=\,\sum_{m_{1}m_{2}t_{1}t_{2}}\,(-)^{j_{1}-m_{1}+1/2-t_{1}}  
 \left(\begin{array}{ccc}
 j_{1} & K & j_{2}\\
 -m_{1} & \kappa & m_{2}\end{array}\right)\left(\begin{array}{ccc}
 1/2 & S & 1/2\\
 -t_{1} & \sigma & t_{2}\end{array}\right)\, a_{1}^{\dagger}a_{2}
 \f] 
In return we use left-to-right order 
\f[
 ({\cal M}_{K\,\kappa,\, S\,\sigma}^{({\bf 1}\,{\bf 2})})^{\dagger}=\sum_{12} q_{21} a_{1}^{\dagger}a_{2} \quad {\cal M}_{K\,\kappa,\, S\,\sigma}^{({\bf 1}\,{\bf 2})}=\sum_{12} q_{12} a_{1}^{\dagger}a_{2}
\f] STensor \f$ \{2 1 \}\rightarrow q_{21} \f$ is returned. No forward ordering, since operator is not generally Hermitian. When operator is used, creation is assumed.
*/ 
int AddMultipole(
  tav::STensor<spsint,double> &ESP,  ///< Output map of e[12] elements
  SingleParticleStates & sps, ///< s.p. states matrix
  const int *ijk ///< list 1 2 LM TS
              ) 
  {

    double htmp;
    bool phase;
    spsint a[2]; 
    spsint b[2]; //whitch order
//order of labels 0123 but operators 
//a0+ a1 a3+ a2
//here we implement a brute force scan
    for (a[0]=0;a[0]<sps.size();a[0]++) {
      if (sps[a[0]].level!=ijk[0]) continue;
      for (a[1]=0;a[1]<sps.size();a[1]++) {
        if (sps[a[1]].level!=ijk[1]) continue;
//spin 
        //enforce projection
            if (sps[a[0]].Jz-sps[a[1]].Jz!=2*ijk[3]) continue;
            if (sps[a[0]].Tz-sps[a[1]].Tz!=2*ijk[5]) continue;
        
            htmp=	tav::ThreeJSymbol(sps[a[0]].J, -sps[a[0]].Jz, 
                                          2*(ijk[2]), (sps[a[0]].Jz-sps[a[1]].Jz),
                                          sps[a[1]].J, sps[a[1]].Jz);
//isospin
            htmp*=	tav::ThreeJSymbol(sps[a[0]].T, -sps[a[0]].Tz, 
                                          2*(ijk[4]), (sps[a[0]].Tz-sps[a[1]].Tz),
                                          sps[a[1]].T, sps[a[1]].Tz);

            if (tav::Abs(htmp)<1E-9) continue;

// cout<<int(a[0])<<" htmp="<<htmp<<" : "<<int(a[0])<<" "<<int(a[1])<<" "<<int(a[2])<<" "<<int(a[3])<<endl;
//	Pause();
//phase 

            phase=
                (((sps[a[0]].J-sps[a[0]].Jz + sps[a[0]].T-sps[a[0]].Tz )>>1)&1);   


 
            b[0]=a[1];b[1]=a[0];  
//enter data 
            ESP[b]+=(phase ? -htmp : htmp); 
   
          }} 
          return 1;
  } //end of function 
  
int PPInteraction2sp(
tav::STensor<spsint,double> &VSP,  ///< Output map of elements V[1234]  (12)
SM::SingleParticleStates & sps, ///< s.p. states matrix
tav::STensor<int,double> &JV, ///< Input spherical interaction
double scaling=1.0   ///<scaling of the matrix element [1.0]
                    )
{
int flag=0;
  tav::STensor<int,double>::iterator iit;
  //check if we are using Tz by looking for negative
  for (iit=JV.begin(); iit!=JV.end();iit++) if ((iit->first[5])<0) {flag=1; break;}
  //no negatives use normal
  if (flag==0) {
  for (iit=JV.begin(); iit!=JV.end();iit++)  
    AddPPVsp(VSP, sps, iit->first, scaling*(iit->second));}
	//use Tz formalism
  else {
  cerr<<"Isospin not conserved "<<endl;
  for (iit=JV.begin(); iit!=JV.end();iit++)  {
    AddPPTzVsp(VSP, sps, iit->first, scaling*(iit->second));
	//cerr<<" Incert matrix element "<< scaling*(iit->second)<<endl;
	}
	}
return 1;
}

/*!\breif transfer interaction in PN format, used in JISP format
Use isospin but cut by pair projection
*/
int PPInteractionPN2sp(
tav::STensor<spsint,double> &VSP,  ///< Output map of elements V[1234]  (12)
SM::SingleParticleStates & sps, ///< s.p. states matrix
tav::STensor<int,double> &JV, ///< Input spherical interaction
int pnstate ///<extra condition on pn-state 1-NN, 0-PN, -1 PP
                    )
{
  tav::STensor<int,double>::iterator iit;
  if (std::abs(pnstate)<2) {
  for (iit=JV.begin(); iit!=JV.end();iit++)  
    AddPPVsp(VSP, sps, iit->first, (iit->second),pnstate);}
    else {
    for (iit=JV.begin(); iit!=JV.end();iit++)  
    AddPPVsp(VSP, sps, iit->first, (iit->second));
    }
    
return 0;
} 


} //end of SM
using SM::SymmetryMBOperator;

#endif

