/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


#ifndef __HO_cxx__
#define __HO_cxx__

/*!
\author Alexander Volya  <http://www.volya.net>
\file HO.cxx
\brief Harmonic Oscillator Wave-functions
\par Definitions
n-main quantum number (number of nodes)
function HORadialIntegral calculates r^L expectation vlaue in HO radial
funtions 
*/

#include <cmath>
#include <tav/basic.cxx>
//#include <mav/gamma.cxx>
//#include <debug.h>
namespace SM{


///Factorial in integer
double SiFactorial(int n) {
  //this in normal factorial
  double ret=1.0;
  for (int i=1;i<=n;i++) ret*=i;
  return ret;}
///Factorial in double
double DiFactorial(int n) {
  //this is double factorial odd only, even gives zero
  double ret=1.0;
  for (int i=n%2;i<=n;i+=2) ret*=i;
  return ret;}
double LogSiFactorial(int n) {
  //this in normal factorial
  double ret=0.0;
  for (int i=1;i<=n;i++) ret+=log(i);
  return ret;}
///Factorial in double
double LogDiFactorial(int n) {
  //this is double factorial odd only, even gives zero
  double ret=0.0;
  for (int i=n%2;i<=n;i+=2) ret+=log(i);
  return ret;}
    /*!
\brief Harmonic oscillator radial integral \f$ \langle n_1 l_1|r^\lambda|n_2 l_2 \rangle \f$

     This is a radial integral for HO wave functions from Lawson book.      
     It is only valid for (l_1 + l_2 + L) even.
     See equation 1.11a page 7 in Lawson
     
     Phase convention is set by bool variable: true (default).
        False=positive at infinity, result is always positive. Phys. Rev. C70, 044005 (2004)
        (the same as in Phys. Lett. B644, 33 (2007)),
        i.e. with extra (-1)^n.
     
        True= wf is positive near origin (Lawson book, page 7). I.Talmi phase convention for oscillator functions;
        (Square root, exponent, laguerre polynomial)
        Tests have shown different definitions could be used
     the Log-based version is overflow protected
     */
    double HORadialIntegral(int L, ///< power \f$ \lambda \f$
                            int n1, ///< first n
                            int l1, ///< first l
                            int n2, ///< second n
                            int l2, ///< second l
                            bool PhaseConvention=true ///< extra phase factor
    
    ) {
        double sterm=1;
        //Determine range of q
        int qmin;
        int qmax;
        qmin=tav::Max(0,(n1-(l2-l1+L)/2), (n2-(l1-l2+L)/2));
        qmax=tav::Min(n1,n2); //we do not devide by zero for qmax
        //  cout<<"qmin= "<<qmin<<endl;
        // cout<<"qmax= "<<qmax<<endl;
        if (qmin>qmax) return 0;
        double numerator=LogDiFactorial(l1+l2+L+2*qmin+1); //this is first term
        double denominator=qmin*log(2.0)+LogSiFactorial(qmin)+LogSiFactorial(n1-qmin)+LogSiFactorial(n2-qmin)+
        LogSiFactorial(qmin+(l1-l2+L)/2-n2)+LogSiFactorial(qmin+(l2-l1+L)/2-n1); //this is term
        double sum=exp(numerator-denominator);
        //  cout<<qmin<<" "<<sum<<endl;
        // numerator=denominator=1;
        for (int q=qmin+1;q<=qmax;q++) {
            numerator+=log(l1+l2+L+2*q+1)+log(n1-q+1)+log(n2-q+1); //decreasing factorial
            denominator+=log(2.*q)+log(q+(l1-l2+L)/2-n2)+log(q+(l2-l1+L)/2-n1);
            //if (q!=n1) denominator*=(n1-q);
            //if (q!=n2) denominator*=(n2-q); //this is last one 0!=1
            //  cout<<q<<" "<<numerator/denominator<<" "<<numerator<<" "<<denominator<<endl;
            sum+=exp(numerator-denominator);
        }
        sum=log(sum);
        //now the sum is debugged, finalize
        sum-=0.5*LogDiFactorial(2*n1+2*l1+1);
        sum-=0.5*LogDiFactorial(2*n2+2*l2+1);
        sum+=0.5*(LogSiFactorial(n1)+LogSiFactorial(n2));
        sum+=((n1+n2-L)/2.0)*log(2.0);
        sum+=LogSiFactorial((l2-l1+L)/2);
        sum+=LogSiFactorial((l1-l2+L)/2);
        //the default (true) positive at origin phase convention has extra (-1)^{n-n'}
        //if ( PhaseConvention && ( (n1+n2) % 2 != 0 )) sum = -sum;
        if ( PhaseConvention && ( (n1+n2)&1 )) return -exp(sum);
        else return exp(sum);
    }
    double HORadialIntegral1(int L, ///< power \f$ \lambda \f$
                            int n1, ///< first n
                            int l1, ///< first l
                            int n2, ///< second n
                            int l2, ///< second l
                            bool PhaseConvention=true ///< extra phase factor
    
    ) {
        double sterm=1;
        //Determine range of q
        int qmin;
        int qmax;
        qmin=tav::Max(0,(n1-(l2-l1+L)/2), (n2-(l1-l2+L)/2));
        qmax=tav::Min(n1,n2); //we do not divide by zero for qmax
        //  cout<<"qmin= "<<qmin<<endl;
        // cout<<"qmax= "<<qmax<<endl;
        if (qmin>qmax) return 0;
        double numerator=DiFactorial(l1+l2+L+2*qmin+1); //this is first term
        double denominator=pow(2., qmin)*SiFactorial(qmin)*SiFactorial(n1-qmin)*SiFactorial(n2-qmin)*
        SiFactorial(qmin+(l1-l2+L)/2-n2)*SiFactorial(qmin+(l2-l1+L)/2-n1); //this is term
        double sum=numerator/denominator;
        //  cout<<qmin<<" "<<sum<<endl;
        // numerator=denominator=1;
        for (int q=qmin+1;q<=qmax;q++) {
            numerator*=(l1+l2+L+2*q+1)*(n1-q+1)*(n2-q+1); //decreasing factorial
            denominator*=2.*q*(q+(l1-l2+L)/2-n2)*(q+(l2-l1+L)/2-n1);
            //if (q!=n1) denominator*=(n1-q);
            //if (q!=n2) denominator*=(n2-q); //this is last one 0!=1
            //  cout<<q<<" "<<numerator/denominator<<" "<<numerator<<" "<<denominator<<endl;
            sum+=numerator/denominator;
        }
        //now the sum is debugged, finalize 
        sum/=sqrt(DiFactorial(2*n1+2*l1+1));
        sum/=sqrt(DiFactorial(2*n2+2*l2+1));
        sum*=sqrt(SiFactorial(n1)*SiFactorial(n2));
        sum*=pow(2.0, (n1+n2-L)/2.0);
        sum*=SiFactorial((l2-l1+L)/2);
        sum*=SiFactorial((l1-l2+L)/2);
        //the default (true) positive at origin phase convention has extra (-1)^{n-n'}
        //if ( PhaseConvention && ( (n1+n2) % 2 != 0 )) sum = -sum;
        if ( PhaseConvention && ( (n1+n2)&1 )) return -sum;
        else return sum;
    }

#ifdef __XMLPARSER__CXX__

    int HORadialIntegral(tav::Tensor<2,double> &RX, xml::XMLNode xValence, int L) {
   int levels=xValence.nChildNode(xsm::valencespace_orbital);
   if (levels==0) FatalError("There are no s.p. levels");
    RX.FullResize(levels,levels);
	std::vector<int> mqn; //main quantum number
	std::vector<int> lll; //orbital quantum number
    xml::XMLNode xtmp;
	for (int i=0; i<levels; i++) {
	xtmp=xValence.getChildNode(xsm::valencespace_orbital,i);

	int itmp; //temporary integer
	xml::getXMLData(itmp,"n", xtmp); mqn.push_back(itmp);
	xml::getXMLData(itmp,"l", xtmp); lll.push_back(itmp);
	//x[i].tag=SPLabel(x[i].N, x[i].L, x[i].J);
	}
	for (int i=0;i<levels;i++) for (int j=0;j<levels;j++) {
   RX[i][j]=SM::HORadialIntegral(L,mqn[i],lll[i], mqn[j],lll[j]);
 }
return 1;
}
#endif //end xml parser part

    
    /*!
     \brief 
     Calculates an integral of r^k between two harmonic oscillator radial wavefunctions 
     (both with an oscillator parameter ν=mω/(2 hbar)). 
     This does not assume that there is a Y_{km}
     associated with the operator so it works for any L
     phase positive at origin (PhaseConventio=true;)
       both versions are really unstable for large values of parameters
     */
   double r_to_the_k_HO_wavefunctions_log(int n1,int l1,int k,int n2,int l2,double nu=0.5)
     {
         //    if (l1!=l2) continue;
        
         //double NF,NF2;
         //double pi=3.141592653589793238462643;
         
//         NF=sqrt(2*pow(2*nu,double(l1+1.5))*tgamma(n1+1)/tgamma(n1+l1+1.5));
//         NF2=sqrt(2*pow(2*nu,double(l2+1.5))*tgamma(n2+1)/tgamma(n2+l2+1.5));
         double lNF=0.5*(log(2.)+double(l1+1.5)*log(2.*nu)+lgamma(n1+1)-lgamma(n1+l1+1.5));
         double lNF2=0.5*(log(2.)+double(l2+1.5)*log(2.*nu)+lgamma(n2+1)-lgamma(n2+l2+1.5));
         double lden=log(2.0)+((l1+l2+k+3)/2.)*log(2.*nu);
         //double prod2=1;
         double sum=0;
         for (int p=0;p<=n2;p++)
         {
             
             
             double sum2=0;
             for (int q=0;q<=n1;q++)
             {
                 double lprod2=lgamma(n1+l1+1.5)-(lgamma(q+l1+1.5));
                 lprod2-=lgamma(q+1)+lgamma(n1-q+1)-(lgamma(p+q+(l1+l2)/2.+(k+3)/2.));
                 if ((q%2)==0) sum2+=exp(lprod2); else sum2-=exp(lprod2);
//                 double prod2=tgamma(n1+l1+1.5)/double(tgamma(q+l1+1.5));
//                 prod2/=tgamma(q+1)*tgamma(n1-q+1)/double(tgamma(p+q+(l1+l2)/2.+(k+3)/2.));
//                 if ((q%2)==0) sum2+=prod2; else sum2-=prod2;
             }
             //double prod=1;
             //if ((p%2)==0) prod=1;
             //double prod=1.0;
             double lprod=lgamma(n2+l2+1.5)-(lgamma(p+l2+1.5));
             lprod-=lgamma(p+1)+lgamma(n2-p+1);
//             double prod=1.0;
//             prod*=tgamma(n2+l2+1.5)/(tgamma(p+l2+1.5));
//             prod/=tgamma(p+1)*tgamma(n2-p+1);
             if ((p%2)==0) sum+=exp(lprod+lNF+lNF2-lden)*sum2; else sum-=exp(lprod+lNF+lNF2-lden)*sum2;
         }
        // sum/=2*pow(2*nu,(l1+l2+k+3)/2.);
         //sum*=exp(lNF+lNF2-lden);
         return sum;
         
     }
//original non-log functions
    double r_to_the_k_HO_wavefunctions(int n1,int l1,int k,int n2,int l2,double nu=0.5)
    {
        //    if (l1!=l2) continue;
        double sum=0;
        double prod=1.;
        double NF,NF2;
        //double pi=3.141592653589793238462643;
        
        NF=sqrt(2*pow(2*nu,double(l1+1.5))*tgamma(n1+1)/tgamma(n1+l1+1.5));
        NF2=sqrt(2*pow(2*nu,double(l2+1.5))*tgamma(n2+1)/tgamma(n2+l2+1.5));
        double sum2=0,prod2=1;
        
        for (int p=0;p<=n2;p++)
        {
            
            prod=-1;
            sum2=0;
            for (int q=0;q<=n1;q++)
            {
                prod2=-1;
                if ((q%2)==0) prod2=1;
                prod2*=tgamma(n1+l1+1.5)/double(tgamma(q+l1+1.5));
                prod2/=tgamma(q+1)*tgamma(n1-q+1)/double(tgamma(p+q+(l1+l2)/2.+(k+3)/2.));
                sum2+=prod2;
            }
            if ((p%2)==0) prod=1;
            prod*=sum2;
            prod*=tgamma(n2+l2+1.5)/(tgamma(p+l2+1.5));
            prod/=tgamma(p+1)*tgamma(n2-p+1);
            sum+=prod;
        }
        sum/=2*pow(2*nu,(l1+l2+k+3)/2.);
        sum*=(NF*NF2);
        return sum;
        
    }




/*
 Function used in evaluation of radial HO oscillator integrals
 T. A. Brody, G. Jacob, M. Moshinsky, Matrix elements in nuclear shell theory, Nuclear Physics 17 (1960) 16--29. Eq. 23,
 see also discussion
 M. Moshinsky, Y. F. Smirnov, The Harmonic Oscillator In Modern Physics, (1996) 2.5
 tried long double and not very good
 */
//typedef long double double;
double HO_B(
int n1, ///< n1
int l1, ///<l1
int n2, ///<n2
int l2, ///< l2
int p  ///< parameter p
){
    if ((l1+l2)&1) return 0; //the l1 and l2 must be of the same pairity
    //limits on p from eq 19 (l1+l2)/2<=p<=n1+n2+(l1+l2)/2 but do not inforce them here
   // if (p<(l1+l2)/2) return 0;
   // if (p>n1+n2+(l1+l2)/2) return 0;
    double logprefator=lgamma(2*p+2)-lgamma(p+1)-(n1+n2)*log(2);
    //cout<<"px="<<exp(logprefator)<<endl;
    logprefator+=0.5*(lgamma(n1+1)+lgamma(n2+1)+lgamma(2*n1+2*l1+2)+
                      lgamma(2*n2+2*l2+2)-lgamma(n1+l1+1)-lgamma(n2+l2+1));

    int kmin;
    if (p-(l1+l2)/2-n2>0) kmin= p-(l1+l2)/2-n2; else kmin=0;  ///eq, 21
    int kmax;
    if (p-(l1+l2)/2>n1) kmax=n1; else kmax=p-(l1+l2)/2;
    //cout<<"kminmax="<<kmin<<" "<<kmax<<endl;
        double logterm=lgamma(l1+kmin+1)+lgamma(p-(l1-l2)/2-kmin+1);
        logterm-=lgamma(kmin+1)+lgamma(2*l1+2+2*kmin)+lgamma(n1-kmin+1);
        logterm-=lgamma(2*p-l1+l2+2-2*kmin)+lgamma(n2-p+(l1+l2)/2+kmin+1);
        logterm-=lgamma(p-(l1+l2)/2-kmin+1);
        //add prefactor
        logterm+=logprefator;
        double term=exp(logterm); //this is a first term in the sum including prefactor
        double sum=term;
        //double testterm=tgamma(p-(l1-l2)/2-kmin+1);
        //double testterm=tgamma(2*p-l1+l2+2-2*kmin);
    	//cout<<kmin<<" "<<p-(l1+l2)/2-kmin<<" "<<term<<" tt0="<<testterm<<endl;
        for (int k=kmin+1;k<=kmax;++k) {
            term*=double(l1+k)/double(k);
            term/=double(p-(l1-l2)/2-k+1);
            term/=double(2*l1+1+2*k)*double(2*l1+1+2*k-1);
            term*=double(n1-k+1);
            term*=double(2*p-l1+l2+2-2*k)*double(2*p-l1+l2+3-2*k);
            term/=double(n2-p+(l1+l2)/2+k);
            term*=double(p-(l1+l2)/2-k+1);
           // testterm/=double(p-(l1-l2)/2-k+1);
            //testterm/=double(2*p-l1+l2+1-2*k+1)*double(2*p-l1+l2+1-2*k+2);
            //cout<<k<<" "<<p-(l1+l2)/2-k<<" "<<term<<" tt="<<testterm<<" "<<tgamma(2*p-l1+l2+2-2*k)<<endl;
            sum+=term;
        }
    //cout<<"sum="<<sum<<" prefector="<<exp(logprefator)<<endl;
    //sum*=exp(logprefator);
    if ((p-(l1+l2)/2)&1) return -sum;  else return sum;
}//end function HO_B
//simplified function, direct gammas for factorials, no recursive division
double HO_B1(
int n1, ///< n1
int l1, ///<l1
int n2, ///<n2
int l2, ///< l2
int p  ///< parameter p
){
    if ((l1+l2)&1) return 0; //the l1 and l2 must be of the same pairity
    //limits on p from eq 19 (l1+l2)/2<=p<=n1+n2+(l1+l2)/2 but do not inforce them here
   // if (p<(l1+l2)/2) return 0;
   // if (p>n1+n2+(l1+l2)/2) return 0;
    double logprefator=lgamma(2*p+2)-lgamma(p+1);
    double logprefator1=-(n1+n2)*log(2);
    logprefator1+=0.5*(lgamma(n1+1)+lgamma(n2+1)+lgamma(2*n1+2*l1+2)+
                      lgamma(2*n2+2*l2+2)-lgamma(n1+l1+1)-lgamma(n2+l2+1));
    int kmin;
    if (p-(l1+l2)/2-n2>0) kmin= p-(l1+l2)/2-n2; else kmin=0;  ///eq, 21
    int kmax;
    if (p-(l1+l2)/2>n1) kmax=n1; else kmax=p-(l1+l2)/2;
    //cout<<"kminmax="<<kmin<<" "<<kmax<<endl;
    double sum=0.0;
    for (int k=kmin;k<=kmax;++k) {
        double logterm=lgamma(l1+k+1)+lgamma(p-(l1-l2)/2-k+1);
        logterm-=lgamma(k+1)+lgamma(2*l1+2+2*k)+lgamma(n1-k+1);
        logterm-=lgamma(2*p-l1+l2+2-2*k)+lgamma(n2-p+(l1+l2)/2+k+1);
        logterm-=lgamma(p-(l1+l2)/2-k+1);
        //add prefactor
        logterm+=logprefator;
       // cout<<"k="<<k<<" "<<exp(logterm)<<endl;
        sum+=exp(logterm);
    }
    if ((p-(l1+l2)/2)&1) return -sum*exp(logprefator1);  else return sum*exp(logprefator1);
}//end function HO_B1
/*equation 19 with coefficients from
2.9 in M. Moshinsky, Y. F. Smirnov, The Harmonic Oscillator In Modern Physics, (1996) 2.5
 after n=15 precision becomes unpredictable
 */
double HORadialIntegralPower(int L, ///< power \f$ \lambda \f$
  int n1, ///< first n
  int l1, ///< first l
  int n2, ///< second n
  int l2 ///< second l
) {
    double sum=0.0;
    for (int p=(l1+l2)/2;p<=n1+n2+(l1+l2)/2;++p) { sum+=HO_B(n1,l1,n2,l2,p)/tgamma(p+1.5)*tgamma(p+1.5+0.5*L);
        //cout<<p<<" "<<tgamma(p+1.5+0.5*L)/tgamma(p+1.5)<<" "<<HO_B1(n1,l1,n2,l2,p)<<endl;
    }
        return sum;
}//end function HORadialIntegralPower

//this program computes a power r^L but changing the order of sums. Still not good convergence.
double HORadialIntegralPower1(int L, ///< power \f$ \lambda \f$
  int n1, ///< first n
  int l1, ///< first l
  int n2, ///< second n
  int l2 ///< second l
) {
    double sum=0.0;
    
     if ((l1+l2)&1) return 0; //the l1 and l2 must be of the same pairity
     //limits on p from eq 19 (l1+l2)/2<=p<=n1+n2+(l1+l2)/2 but do not inforce them here
    // if (p<(l1+l2)/2) return 0;
    // if (p>n1+n2+(l1+l2)/2) return 0;
    double logexternal=-(n1+n2)*log(2)+0.5*(lgamma(n1+1)+lgamma(n2+1)+lgamma(2*n1+2*l1+2)+
    lgamma(2*n2+2*l2+2)-lgamma(n1+l1+1)-lgamma(n2+l2+1));
    for (int k=0;k<=n1;++k) {
    double sumk=0;
    double logk=lgamma(l1+k+1)-lgamma(k+1)-lgamma(2*l1+2+2*k)-lgamma(n1-k+1);
    for (int p=(l1+l2)/2;p<=n1+n2+(l1+l2)/2;++p) {
 
                  if (n2-p+(l1+l2)/2+k<0) continue; //skip until we reach correct kmin
                 if (p-(l1+l2)/2-k<0) continue;
        double logprefator=lgamma(2*p+2)-lgamma(p+1)+lgamma(p+1.5+0.5*L)-lgamma(p+1.5);//
        
         double logterm=lgamma(p-(l1-l2)/2-k+1);
         logterm-=lgamma(2*p-l1+l2+2-2*k)+lgamma(n2-p+(l1+l2)/2+k+1);
         logterm-=lgamma(p-(l1+l2)/2-k+1);
         //add prefactor
         logterm+=logprefator;
        // cout<<k<<" "<<p-(l1+l2)/2-k<<" "<<exp(logterm)<<endl;
         //sum+=exp(logterm);
       // cout<<"exp "<<exp(logterm)<<endl;
         if ((p-(l1+l2)/2)&1) sumk-=exp(logterm);  else sumk+=exp(logterm);
         //if (k==p-(l1+l2)/2) break; //stop if we reach zero on next step
     }
       sum+=sumk*exp(logk);
    }
        return sum*exp(logexternal);
}//end function HORadialIntegralPower

} //end namespace

#endif
