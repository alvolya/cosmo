/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


#ifndef __R__TO__THE__K__INTERACTION__CXX__
#define __R__TO__THE__K__INTERACTION__CXX__

#include <mav/moshinsky.cxx>
#include <mav/NineJSymbol.cxx>
#include <mav/ThreeJSymbol.cxx>
//#include <kmav/r_to_the_k.cxx>
using namespace mav;
#include <sm/HO.cxx>

#include <cmath>
#include <iostream>
using namespace std;

/*! \brief
 Calculates matrix elements of  interaction that is |r_i - r_j|^k
 */
double rToTheK_IT(int n1,int l1,double j1,int n2,int l2,double j2,int n3, int l3,double j3,int n4,int l4, double j4,double J,int T,int input_k)
{
    //    cout<<T<<" ";
    double rtrn=1.;
    if (n1==n2 && l1==l2 && j1==j2)
        rtrn/=sqrt(2.);
    if (n3==n4 && l3==l4 && j3==j4)
        rtrn/=sqrt(2.);
    rtrn*=sqrt( (2*j1+1) * (2*j2+1) * (2*j3+1) * (2*j4+1));
    int qi=2*(n1+n2)+l1+l2;
    int qj=2*(n3+n4)+l3+l4;
    //    if (qi!=qj) return 0;
    //    int qmin=( (qi<qj) ? (qi) : (qj) );
    int nphase=1;
    
    int lambdamin=( (abs(l1-l2)<abs(l3-l4)) ? (abs(l3-l4)) : (abs(l1-l2)) );
    int lambdamax=( (abs(l1+l2)>abs(l3+l4)) ? (abs(l3+l4)) : (abs(l1+l2)) );
    double sum=0,prod=0;
    double sq2=pow(M_SQRT2,input_k);//This is here because |r_1 - r_2|^k =Sqrt(2)^k * r^k so it saves some time.
    for (int lambda=lambdamin;lambda<=lambdamax;lambda++)
    {
        for (int S=0;S<=1;S++)
        {
            for (int n=0;n<=qi/2;n++)
            {
                for (int l=0;l<=qi;l++)
                {
                    if ((l+S+T)%2==0) continue;
                    
                    for (int N=0;N<=qi/2;N++)
                    {
                        for (int L=0;L<=qi;L++)
                        {
                            if (2*n+2*N+L+l != qi) continue;
                            
                            prod=2*(2*lambda+1)*(2*S+1);
                            if (prod==0) continue;
                            //                            cout<<prod<<" > NJ1: ";
                            prod*=NineJSymbol(l1,l2,lambda,0.5,0.5,S,j1,j2,J);
                            if (std::abs(prod)<1.e-9) continue;
                            //                            cout<<prod<<" > NJ2: ";
                            prod*=NineJSymbol(l3,l4,lambda,0.5,0.5,S,j3,j4,J);
                            if (std::abs(prod)<1.e-9) continue;
                            //                            cout<<prod<<" ";
                            prod*=MoshinskyBrackets(n,l,N,L,n1,l1,n2,l2,lambda,1.);
                            if (std::abs(prod)<1.e-9) continue;
                            int np=(qj-l-L)/2 -N;
                            prod*=MoshinskyBrackets(np,l,N,L,n3,l3,n4,l4,lambda,1.);
                            if (std::abs(prod)<1.e-9) continue;
                            //                            prod*=Integrate(n,l,np,l,sign,nu);
//                            if ((n%2)!=(np%2)) nphase=-1;//This is here  because HO is defined with (-)^n phase wf's
//                            else nphase=1;
//                            prod*=SM::HORadialIntegral(input_k,n,l,np,l)*nphase;//over sqrt(2) because r=(r1-r2)/sqrt2 to have moshinsky coefficients
                            prod*=SM::r_to_the_k_HO_wavefunctions(n,l,input_k,np,l,1./2.);
                            //                            nphase=1;
                            sum+=(prod);
                        }
                    }
                }
            }
        }
    }
    //    cout<<rtrn*sum<<endl;
    return (rtrn*sum)*sq2;
}
//template<class Space>
int GeneraterToTheKINT(tav::STensor<int,double> &VSP, SM::ValenceSpace &x,int k_)
{
    GenerateEmptyINT(VSP,x);//First make empty INT tensor.
    tav::STensor<int,double>::iterator it=VSP.begin();
    int s1,s2,s3,s4,J,T;
    for (;it!=VSP.end();it++)
    {
        s1=it->first[0];s2=it->first[1];
        s3=it->first[2];s4=it->first[3];
        J=it->first[4];T=it->first[5];
        it->second=rToTheK_IT(x[s1].N,x[s1].L,int(x[s1].J)/2.,
                              x[s2].N,x[s2].L,int(x[s2].J)/2.,
                              x[s3].N,x[s3].L,int(x[s3].J)/2.,
                              x[s4].N,x[s4].L,int(x[s4].J)/2.,
                              double(J),T,k_);
        
    }
    
    return 0;
}

#endif
