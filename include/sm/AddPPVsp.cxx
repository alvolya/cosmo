/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*
Attempt to make a faster version of the code
*/
#include <tav/stensor.cxx>
#include <tav/ThreeJSymbolM.cxx>
using tav::ThreeJSymbol;
//#include <sm/QN.cxx>
#include <sm/SingleParticleStates.cxx>
namespace SM{
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
New version where CGC list is generated initially
 */
int AddPPVsp(
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
//if (ijk[0]>ijk[1]) cerr<<"levels not ordered "<<endl;
//Here we assume that for a given level the single-particle indexes proceed sequentially although they are broken by isospin in two regions.
spsint a_min1; //set of initial sps
a_min1=0;

while ((a_min1<sps.size())&&(sps[a_min1].level!=ijk[1])) {a_min1+=spsint(sps[a_min1].J+1);} //find first occurence


    for (a[0]=0;a[0]<sps.size();a[0]++) {
    if (sps[a[0]].level!=ijk[0])
    {a[0]+=spsint(sps[a[0]].J); continue; /*add 2J moving toward end and continue */}
    //Let us work on array of CGC here
    int l1=sps[a[0]].J; int m1=sps[a[0]].Jz;
    int l2=sps[a_min1].J;
    int l3=2*(ijk[4]);
    int m2min = ((-l2> -l3-m1) ? -l2 : -l3-m1); //MAX(-l2,-l3-m1), MAX(a,b) ((a>b)?a:b)
    int m2max = ((l2<l3-m1)? l2: l3-m1); //MIN(l2,l3-m1); MIN(a,b) ((a<b)?a:b)
    double *cgcof = new double[((m2max-m2min)/2+1)];
    if (tav::ThreeJSymbolM (l1,l2,l3,m1, m2min,m2max, cgcof)) {
    delete [] cgcof;
    continue; //CGC is unphysical
    //it is hard to see this happening but just in case
    }
    
    for (a[1]=a_min1;a[1]<sps.size();a[1]++) {
    if (sps[a[1]].level!=ijk[1])
    {a[1]+=spsint(sps[a[1]].J); continue; /*add 2J moving toward end and continue */}
    if (((sps[a[0]].J)&1)&&(a[1]==a[0])) continue; //boson fermion
    int m2=sps[a[1]].Jz;
    if (m2<m2min) continue;
    if (m2>m2max) continue;
    double alCGC=cgcof[(m2-m2min)/2]; //this is the only allowed CGC here
    

    for (a[2]=0;a[2]<sps.size();a[2]++) {
//pairs are NOT forward ordered: we can not do that because of mapping to levels
    if (sps[a[2]].level!=ijk[2])
    {a[2]+=spsint(sps[a[2]].J); continue; /*add 2J moving toward end and continue */}
    for (a[3]=0;a[3]<sps.size();a[3]++) {
    if (sps[a[3]].level!=ijk[3])
    {a[3]+=spsint(sps[a[3]].J); continue; /*add 2J moving toward end and continue */}
        if (((sps[a[2]].J)&1)&&(a[2]==a[3])) continue; //fermion and boson
//spin
    if (sps[a[0]].Jz+sps[a[1]].Jz!=sps[a[3]].Jz+sps[a[2]].Jz) continue;
    if (sps[a[0]].Tz+sps[a[1]].Tz!=sps[a[3]].Tz+sps[a[2]].Tz) continue;
        //cout<<int(a[0])<<" "<<int(a[1])<<" "<<int(a[2])<<" "<<int(a[3])<<endl;
        phase=false;
//        htmp=tav::ClebschGordan(sps[a[0]].J, sps[a[0]].Jz,
//                                   sps[a[1]].J, sps[a[1]].Jz,
//                           2*(ijk[4]), (sps[a[0]].Jz+sps[a[1]].Jz));
        htmp=(((l1+m1+m2+3*l2)>>1)&1) ?
    -sqrt(l3+1.0)*alCGC : sqrt(l3+1.0)*alCGC ;

        
// cout<<sps[a[0]].J/2.0 <<" "<< -sps[a[0]].Jz/2.0<<" "<<
//     double(ijk[4])<<" "<<(sps[a[0]].Jz-sps[a[1]].Jz)/2.0<<" "<<
//     sps[a[1]].J/2.0<<" "<<sps[a[1]].Jz/2.0<<endl;
//cout<<int(a[0])<<" htmpxx="<<htmp<<" : "<<int(a[0])<<" "<<int(a[1])<<" "<<int(a[2])<<" "<<int(a[3])<<endl;
htmp*=    tav::ClebschGordan(sps[a[2]].J, sps[a[2]].Jz,
                                   sps[a[3]].J, sps[a[3]].Jz,
                           2*(ijk[4]), (sps[a[2]].Jz+sps[a[3]].Jz));
//cout<<int(a[0])<<" htmpx="<<htmp<<" : "<<int(a[0])<<" "<<int(a[1])<<" "<<int(a[2])<<" "<<int(a[3])<<endl;
//isospin
 htmp*=tav::ClebschGordan(sps[a[0]].T, sps[a[0]].Tz,
                         sps[a[1]].T, sps[a[1]].Tz,
                             2*(ijk[5]), (sps[a[0]].Tz+sps[a[1]].Tz));
 htmp*=    tav::ClebschGordan(sps[a[2]].T, sps[a[2]].Tz,
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
 
    }}}
    
    delete [] cgcof;
    }
    return 0;
} //end of function
}
