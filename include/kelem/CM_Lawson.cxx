/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


#ifndef __CM__LAWSON__CXX__
#define __CM__LAWSON__CXX__

#include <mav/moshinsky.cxx>
#include <mav/gamma.cxx>
#include <mav/NineJSymbol.cxx>
#include <mav/ThreeJSymbol.cxx>
using namespace mav;
#include <sm/SingleParticleStates.cxx>
#include <sm/sphericalint.cxx>
using namespace SM;
#include <cmath>
#include <iostream>
using namespace std;

//int Delta(int a,int b)
//{
//    if (a==b) return 1;
//    else return 0;
//}
double E_IT(int n1,int l1,double j1,int n2,int l2,double j2,int n3, int l3,double j3,int n4,int l4, double j4,double J,int T)
{
//default spin
    double ss=std::abs(j1-l1);
    int s2=int(2.0*ss+0.0001); //max spin integer
    //cout<<"spin "<<ss<<endl;
    double rtrn=1.;
    if (n1==n2 && l1==l2 && j1==j2)
        rtrn/=sqrt(2.);
    if (n3==n4 && l3==l4 && j3==j4)
        rtrn/=sqrt(2.);
    rtrn*=sqrt( (2*j1+1) * (2*j2+1) * (2*j3+1) * (2*j4+1));
    int qi=2*(n1+n2)+l1+l2;
    int qj=2*(n3+n4)+l3+l4;
    if (qi!=qj) return 0;
    int qmin=( (qi<qj) ? (qi) : (qj) );
    
    int lambdamin=( (abs(l1-l2)<abs(l3-l4)) ? (abs(l3-l4)) : (abs(l1-l2)) );
    int lambdamax=( (abs(l1+l2)>abs(l3+l4)) ? (abs(l3+l4)) : (abs(l1+l2)) );
    double sum=0,prod=0;
    for (int lambda=lambdamin;lambda<=lambdamax;lambda++)
    {
        for (int S=0;S<=s2;S++)
        {
            for (int n=0;n<=qmin/2;n++)
            {
                for (int l=0;l<=qmin;l++)
                {
                    if ((l+S+T)%2==0) continue;
                    
                    for (int N=0;N<=qmin/2;N++)
                    {
                        for (int L=0;L<=qmin;L++)
                        {
                            if (2*n+2*N+L+l != qmin) continue;

                            prod=2*(2*lambda+1)*(2*S+1)*(2*N+L-2*n-l);
                            if (prod==0) continue;
//                            cout<<prod<<" > NJ1: ";
                            prod*=NineJSymbol(l1,l2,lambda,ss,ss,S,j1,j2,J);
                            if (prod==0) continue;
//                            cout<<prod<<" > NJ2: ";
                            prod*=NineJSymbol(l3,l4,lambda,ss,ss,S,j3,j4,J);
                            if (prod==0) continue;
//                            cout<<prod<<" ";
                            prod*=MoshinskyBrackets(n,l,N,L,n1,l1,n2,l2,lambda,1.);
                            if (prod==0) continue;
//                            cout<<prod<<" ";
                            prod*=MoshinskyBrackets(n,l,N,L,n3,l3,n4,l4,lambda,1.);
                            if (prod==0) continue;
                            sum+=(prod);
                        }
                    }
                }
            }
        }
    }
    return (rtrn*sum);
}
//template<class Space>//This is templated now so it can accept both SM::ModelSpace and SM::ValenceSpace
int GenerateCMINT(tav::STensor<int,double> &VSP, SM::ValenceSpace &x)
{
    GenerateEmptyINT(VSP,x);//First make empty INT tensor.
    int s1,s2,s3,s4,J,T;

    std::vector<tav::STensor<int,double>::iterator> itvec(VSP.size(),VSP.begin());//Make Array of iterators
    int k=0;
    for(tav::STensor<int,double>::iterator it=VSP.begin();it!=VSP.end();it++,k++)//Fill Array.
        itvec[k]=it;
#pragma omp parallel for
    for (int i=0;i<itvec.size();i++)
    {
        tav::STensor<int,double>::iterator it=itvec[i];
        s1=it->first[0];s2=it->first[1];
        s3=it->first[2];s4=it->first[3];
        J=it->first[4];T=it->first[5];
        it->second=E_IT(x[s1].N,x[s1].L,int(x[s1].J)/2.,
                        x[s2].N,x[s2].L,int(x[s2].J)/2.,
                        x[s3].N,x[s3].L,int(x[s3].J)/2.,
                        x[s4].N,x[s4].L,int(x[s4].J)/2.,
                        double(J),T);
    }
    //This is non parallel version, left here in case someone suspects this doesn't parallelize properly.
//    tav::STensor<int,double>::iterator it=VSP.begin();
//    for (;it!=VSP.end();it++)
//    {
//        s1=it->first[0];s2=it->first[1];
//        s3=it->first[2];s4=it->first[3];
//        J=it->first[4];T=it->first[5];
//        it->second=E_IT(x[s1].N,x[s1].L,int(x[s1].J)/2.,
//                        x[s2].N,x[s2].L,int(x[s2].J)/2.,
//                        x[s3].N,x[s3].L,int(x[s3].J)/2.,
//                        x[s4].N,x[s4].L,int(x[s4].J)/2.,
//                        double(J),T);
//                        
//    }
    
    return 0;
}
#endif
