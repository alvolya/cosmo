/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file moshinsky.cxx
\brief Moshinsky oscillator brackets
\ingroup gp_mav
*/
#ifndef __MOSHINSKY__CXX__
#define __MOSHINSKY__CXX__


#include <mav/ThreeJSymbol.cxx>
#include <mav/NineJSymbol.cxx>
using std::abs;
//#include <mav/gamma.cxx>
//using mav::LogGamma;
#define LogGamma(x) lgammal(x) //the long double allows more precision
typedef long double dnumber;

/*!
 \brief The purpose of this file is to hold the function Moshinsky Brackets that calculates the Moshinsky Bracket for a given set of quantum numbers using the formula by L. Trlifaj ( http://prc.aps.org/pdf/PRC/v5/i5/p1534_1 )
 
    code starts producing innacurate results at about 20 quanta of energy.
    phase convention true (positive at origin)
 */


/** Main function, 2 diffrent cases l=0, and everything else. The second case is just a sum over the possible l=0 coefficients.
 INPUT: n  (number of nodes for relative coordinate
        l  (relative angular momentum)
        N  (number of nodes for Center of Mass
        L  (Center of Mass angular momentum)
        n1 (number of nodes for particle 1)
        l1 (orbital angular momentum for particle 1)
        n2 (number of nodes for particle 2)
        l2 (orbital angular momentum for particle 2)
        lambda (λ is the total angular momentum the 2 particles are coupled to)
        d  (the mass ratio of the 2 particles d=m2/m1)
 
 */

dnumber MoshinskyBrackets(
int n, ///<n  (number of nodes for relative coordinate)
int l, ///< l  (relative angular momentum)
int N, ///<N  (number of nodes for Center of Mass)
int L, ///<L  (Center of Mass angular momentum)
int n1, ///<n1 (number of nodes for particle 1)
int l1, ///< l1 (orbital angular momentum for particle 1)
int n2, ///< n2 (number of nodes for particle 2)
int l2, ///< l2 (orbital angular momentum for particle 2)
int lambda, ///< lambda (λ is the total angular momentum the 2 particles are coupled to)
dnumber d ///<  d  (the mass ratio d=m2/m1)
)
{
    if ((2*n+2*N+L+l)!=(2*n1+l1+2*n2+l2)) return 0;//quick return.
    if (std::abs(L-l)>lambda || (L+l)<lambda)   return 0;
    if (std::abs(l1-l2)>lambda || (l1+l2)<lambda) return 0;
//    if (((l1+l2)&1)^((L+l)&1)) return 0;//this is fixed by energy conservation but ok.

    
    if (l==0)//case of l=0
    {
        dnumber rtrn=1.;
        if (N&1) rtrn=-.5; else rtrn=0.5;
        rtrn*=pow(1.+d,-l1/2.)*pow(d/dnumber(1.+d),l2/2.);
       // rtrn*=uint_type(n,l,N,L,n1,l1,n2,l2);
        rtrn*=sqrt((2*l1+1)*(2*l2+1)/dnumber(2*L+1));
        dnumber rx=LogGamma(n+1)+LogGamma(n1+1)+2*LogGamma(0.5); //last term log(pi)
        rx-=LogGamma(N+1);
        rx+=LogGamma(n2+1);
        rx+=LogGamma(n1+l1+1.5)-LogGamma(N+L+1.5);
        rx+=LogGamma(n2+l2+1.5)-LogGamma(n+1.5);
        rtrn*=exp(-rx);
        rtrn*=mav::ClebschGordan(l1,0,l2,0,L,0);
        dnumber sum=0.0;
        for (int t1=0;t1<=n1;t1++)
        {
            for (int t2=((N+((L-l1-l2)>>1)-t1)<0 ? 0 : N+((L-l1-l2)>>1)-t1); t2<=n2;t2++)
            //implement limit ((t1+t2-N+(l1+l2-L)/2)>=0)  t2>=N+(L-l1-l2)>>1-t1
            {
              //  if ((t1+t2-N+(l1+l2-L)/2)<0) continue; //number of quanta check


                    dnumber rx1=log(d/(1+d))*t2-log(1+d)*t1;
                    rx1+=LogGamma((l1+l2-L)/2+t1+t2+1.)-LogGamma((l1+l2-L)/2+t1+t2-N+1.);
                    rx1+=LogGamma((l1+l2+L)/2+t1+t2+1.5);
                    rx1-=LogGamma(t1+1);
                    rx1-=LogGamma(t2+1);
                    rx1-=LogGamma(n1-t1+1);
                    rx1-=LogGamma(n2-t2+1);
                    rx1-=LogGamma(t1+l1+1.5);
                    rx1-=LogGamma(t2+l2+1.5);
                    if ((t1+t2)&1) sum-=exp(rx1+1.5*rx); else sum+=exp(rx1+1.5*rx); //parity check x&1 true for odd

                
            }
        }
//        cout<<rtrn<<"\t"<<sum<<endl;
        rtrn*=sum;
        return rtrn;
        
    }
    else
    {
        dnumber rtrn;
        if ((n1+n2+n+lambda+N)&1) rtrn=-1.0; else rtrn=1.0;
        rtrn*=(2*l+1);
        dnumber lgrt=0;// part of the log term to be combined in the sum
        lgrt+=-0.5*l*log(1.+d);
        lgrt+=0.5*(LogGamma(n1+1)+LogGamma(n2+1)+log(2*L+1.));
        lgrt+=0.5*(LogGamma(n+1.5)+LogGamma(n1+l1+1.5)+LogGamma(n2+l2+1.5)-LogGamma(n+l+1.5));
        dnumber sum=0;
//        cout<<exp(lgrt*0.5)<<endl;getchar();
        int lmax=L;
        if (l1>lmax) lmax=l1;
        if (l2>lmax) lmax=l2;
        //given p1>=0 we must have 2*n1+l1>=la1
         for (int la1=0;la1<=(l);la1++)
        {           
        int la2=l-la1;
        //la1 and la2 are defined, put limits on p:  mqn: 2*n1-la1+l1>=p1 triangle rule  la1+l1>=p1>=|la1-l1|
        //this means abs(la1-l1)<=p1<=l1+min(n1,la1)
        //due to parity p1+la1+l1 must be even, go in steps of 2
            for (int p1=abs(la1-l1);p1<=l1+((n1<la1)? 2*n1-la1: la1) ;p1+=2)
            {
             for (int p2=abs(la2-l2);p2<=l2+((n2<la2)? 2*n2-la2: la2) ;p2+=2)
                {

 //--------------check all conditions---------------------------
   // if( (n1-(la1+p1-l1)/2.)<0) continue;
  //  if( (n2-(la2+p2-l2)/2.)<0) continue;
    //Check Gamma Functions. not needed because the above inequality is stronger 
   // if( (n1+(-la1+p1+l1)/2.+1.5)<0) continue;
   // if( (n2+(-la2+p2+l2)/2.+1.5)<0) continue;
    //clebsch gordan
    //if (l1>(p1+la1) || l1<std::abs(p1-la1)) continue;
   // if (l2>(p2+la2) || l2<std::abs(p2-la2)) continue;
    if (L>(p2+p1) || L<std::abs(p2-p1)) continue; //p1p2 also satisfy parity with L, this probably is never used?
//-end check conditions 
 //initilize prod
                  
                    dnumber prod;
                    if ((la1+(p1+p2+L)/2)&1) prod=-1.0; else prod=1.0;
                    dnumber lprod=0.5*la1*log(d);
                    lprod+=0.5*(log((2*p1+1)*(2*p2+1))+LogGamma(2.*l+1.)-LogGamma(2.*la1+1.)-LogGamma(2.*(l-la1)+1.));
                    prod*=mav::ClebschGordan(p1,0,la1,0,l1,0);
//                    cout<<"h "<<(n1-(la1+p1-l1)/2+1)<<endl;
                    lprod-=0.5*(LogGamma(n1-(la1+p1-l1)/2+1)+LogGamma(n2-(la2+p2-l2)/2+1));
//                    cout<<"h2"<<endl;
                    lprod-=0.5*(LogGamma(n1+(p1+l1-la1)/2+1.5)+LogGamma(n2+(p2+l2-la2)/2+1.5));
                   // prod*=exp(lprod+2*lgrt/3.);
                    prod*=exp(lprod+lgrt);
                    prod*=mav::ClebschGordan(p2,0,la2,0,l2,0);
                    prod*=MoshinskyBrackets(n,0,N,L,(n1-(la1+p1-l1)/2),p1,(n2-(la2+p2-l2)/2),p2,L,d);
                    prod*=mav::NineJSymbol(p1,p2,L,la1,la2,l,l1,l2,lambda);
                    sum+=prod;
//                    cout<<prod<<" "<<exp(lprod+lgrt/2.)<<"\n";
                }
            }
        }
      //   cerr<<"rtrn ="<<rtrn*exp(lgrt)<<endl;
       // cerr<<"sum ="<<sum<<endl;
//        cout<<rtrn<<" "<<sum<<endl;
        rtrn*=sum;
//        rtrn*=exp(lgrt/3.);
        return rtrn; 
    }
}

#endif //__MOSHINSKY__CXX__

