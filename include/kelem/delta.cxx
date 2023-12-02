/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/



#ifndef __delta__interaction__cxx__
#define __delta__interaction__cxx__

#include <cmath>
#include <iostream>
#include <debug.h>
//#include <mav/SixJSymbol.cxx>
#include <tav/ThreeJSymbolM.cxx>
using namespace std;



/** \brief Delta interaction \f$ V(r_i-r_j) =-4\pi \delta(r_i-r_j)\f$
 \f[
 V_{JT}=-\frac{(-1)^{j_1+j_2+j_3+j_4}}{2L+1}
 \sqrt{\frac{(2j_1+1)(2j_2+1)(2j_3+1)(2j_4+1)}
 {(1+\delta_{j_1j_2}\delta_{l_1 l_2})(1+\delta_{j_3j_4}\delta_{l_3 l_4})}} \times \f]
 \f[
 \left [
 \frac{1+(-1)^{T}}{2}(j_1 1/2, j_2 1/2| L 1) (j_3 1/2, j_4 1/2| L 1)+
 (-1)^{l_2+l_4+j_2-j_4}
 \frac{1-(-1)^{L+T+l_3+l_4}}{2}(j_1 1/2, j_2 -1/2| L 0) (j_3 1/2, j_4 -1/2| L 0)
 \right ]
 \f]
 The check for parity is also required by l.
 The result is to be multiplied by a radial integral
 \f[
 \overline{R}=\int R_1(r)R_2(r)R_3(r)R_4(r) r^2 dr
 \f]
 See Lawson, Theory of the Nuclear Shell Model, page 438
 For SDI \f$\overline{R}=(-1)^{n_1+n_2+n_3+n_4} \f$ where n=0,1 is number of nodes.
 
 */
double VDelta(
              int j1, int l1,
              int j2, int l2,
              int j3, int l3,
              int j4, int l4,
              int L, int T
              )
{
    //parity check
    if (((l1+l2+l3+l4)>>1)&1) return 0.0;
    double X;
    
    if (((T)>>1)&1) X=0;
    else X=tav::ClebschGordan(j1,1,j2,1,L,2)*tav::ClebschGordan(j3,1,j4,1,L,2);
    
    double Y;
    if (((T+L+l3+l4)>>1)&1) {
        Y=tav::ClebschGordan(j1,1,j2,-1,L,0)*tav::ClebschGordan(j3,1,j4,-1,L,0);
        if (((l2+l4+j2-j4)>>1)&1) Y=-Y;
    }
    else Y=0;
    double V;
    
    V=(j1+1.)*(j2+1.)*(j3+1.)*(j4+1.);
    if ((j1==j2)&&(l1==l2)) V/=2.0;
    if ((j3==j4)&&(l3==l4)) V/=2.0;
    V=sqrt(V)*(X+Y)/(L+1.);
    if(((j1+j2+j3+j4)>>1)&1) return V; else return -V;
}


//template<class Space>//This is templated now so it can accept both SM::ModelSpace and SM::ValenceSpace
int GenerateSDIINT(tav::STensor<int,double> &VSP, SM::ValenceSpace &x)
{
    GenerateEmptyINT(VSP,x);//First make empty INT tensor.
    tav::STensor<int,double>::iterator it=VSP.begin();
    int s1,s2,s3,s4,J,T;
    for (;it!=VSP.end();it++)
    {
        s1=it->first[0];s2=it->first[1];
        s3=it->first[2];s4=it->first[3];
        J=it->first[4];T=it->first[5];
        it->second = VDelta(  x[s1].J,2*x[s1].L,
                              x[s2].J,2*x[s2].L,
                              x[s3].J,2*x[s3].L,
                              x[s4].J,2*x[s4].L,
                              2*J,2*T);
        if((x[s1].N+x[s2].N+x[s3].N+x[s4].N)&1) it->second=-it->second;
        
    }
    
    return 0;
}

//This is old version. Replaced with GenerateEmptyINT 
///** \brief Generate Surface-Delta Interaction
// */
//int GenerateSDIINT(tav::STensor<int,double> &VSP, SM::ValenceSpace &x) {
//    int s1,s2,s3,s4;
//    int *aa=new int [6];
//    
//    if (x.size()==0)
//    {cerr<<"Error in GenerateSDIINT, you must define space first"<<endl;
//        return 0;}
//    
//    //we first need to count how many pair states we can have
//    cerr<<x.size()<<endl;
//    int tbme=0;
//    for (s2=0;s2<x.size();s2++)
//        for (s1=s2;s1<x.size();s1++) //this will list full pair (s1>= s2) -oxbash
//            for (s4=s2;s4<x.size();s4++) {
//                if (s4==s2) s3=Max(s4,s1); else s3=s4;
//                //this is where the trick is we compare pair by first number s2 s4
//                //but if those are equal we compare using s1 and s3
//                //rmemember s1>=s2 same s3>=s4
//                //from independence of pairs
//                //s4>=s2 and if s4==s2 then s3>=s1
//                for (;s3<x.size();s3++){
//                    //  cerr<<"("<<s1<<" "<<s2<<" )  ("<<s3<<" "<<s4<<" )";
//                    //now we have all pairs and we need to decide what L and T are alowed
//                    for (int L=abs((x[s1].J-x[s2].J)/2);L<=((x[s1].J+x[s2].J)/2);L++)
//                        for (int T=abs((x[s1].T-x[s2].T)/2);T<=((x[s1].T+x[s2].T)/2);T++)
//                        {
//                            //now we need to check if L and t are alowed for both pairs
//                            if ((2*L)>(x[s3].J+x[s4].J)) continue; //bad triangle rule
//                            if ((2*L)<abs(x[s3].J-x[s4].J)) continue; //bad triangle rule
//                            if ((2*T)>(x[s3].T+x[s4].T)) continue; //bad triangle rule
//                            if ((2*T)<abs(x[s3].T-x[s4].T)) continue; //bad triangle rule
//                            //now if we are on one level I need to check parity
//                            if ((s1==s2)&&(((x[s1].J-L+x[s1].T-T)%2)==0)) continue;
//                            if ((s3==s4)&&(((x[s3].J-L+x[s3].T-T)%2)==0)) continue;
//                            if ((x[s1].L+x[s2].L+x[s3].L+x[s4].L)%2==1) continue; //parity
//                            tbme++;
//                            // cerr<<"("<<s1<<" "<<s2<<" )  ("<<s3<<" "<<s4<<" )"<<" "<<L<<" "<<T<<endl;
//                            
//                            aa[0]=s1; aa[1]=s2; aa[2]=s3; aa[3]=s4; aa[4]=L; aa[5]=T;
//                            double htmp=
//                            VDelta(x[s1].J,2*x[s1].L,x[s2].J,2*x[s2].L,x[s3].J,2*x[s3].L,x[s4].J,2*x[s4].L,2*L, 2*T);
//                            if (fabs(htmp)<1E-10) continue;
//                            if((x[s1].N+x[s2].N+x[s3].N+x[s4].N)&1) VSP[aa]=-htmp; else VSP[aa]=htmp;
//                            
//                        }//this is end cycle over T & L
//                }}
//    //this goes over pair (s3>=s4)
//    /*let us define pair comparison so not to double count pairs
//     (s1,s2) (s3,s4)  if s2>s4 then P12>P34
//     
//     */
//    // cerr<<tbme<<endl;
//    // DBGPrint(tbme);
//    delete [] aa;
//    return tbme;
//}
#endif

