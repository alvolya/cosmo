/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/** 
\author Alexander Volya  <http://www.volya.net>
\file gamma.cxx
\brief Special functions using double: Gamma, Beta, Harmonic Series, BinomialCoefficient .
*/
#ifndef __GAMMA_CXX__
#define __GAMMA_CXX__

#include <debug.h>
#include <cmath>
using std::atan;
using std::fabs;
#include <constants.cxx>
using phy::alpha;

#ifndef PI
#define PI 3.141592653589793238462643
#endif

namespace mav{
/**
 \namespace mav
 \brief Non-Templated Mathematical Functions
 */
/** Coulomb phase for set of momenta \f$\phi_c=\arg [\Gamma(1+l+i\eta)] \f$
  */
double CoulombPhase(double eta, ///<Coulimb parameter 
                    int l, ///<angular momentum  
                    double *slc ///<array with phase shifts for L=0...l
                   ){
//c *** sigmal(eta,l)= Coulomb phase = arg[Gamma(l+1+i*eta)]
//c *** all Coulomb phase from l=0 to given l are placed in slc(l)
//c *** so that sigma_l=slc(l+1) (in radians)
      double C=0.577215664901533;
      double  eps=1.0E-13;
      double sigma0=-C*eta;
      double rl=double(l);
      double rm=0.0;
      double c1m,cm, sigmal;
   l2:  
      c1m=eta/(1.0+rm);
      cm=c1m-atan(c1m);
      if(fabs(cm)<eps) goto l1;
      sigma0=sigma0+cm;
      rm=rm+1.0;
      goto l2;
   l1: 
      slc[0]=sigma0;
      if(l==0) return sigma0;
      for (int m=1;m<=l;m++)
      slc[m]=slc[m-1]+atan(eta/double(m));
      return slc[l];
}
  
  /** Coulomb phase \f$\phi_c=\arg [\Gamma(1+l+i\eta)] \f$
  
The test code is  columbphase_test.cpp
Coulomb parameter 
  \f[ 
  \eta=e^2 \frac{zZ}{\hbar v}=\alpha \frac {c}{v}
  \f]
  */
double CoulombPhase(double eta, ///<Coulimb parameter \f$ \eta=\alpha/k \f$
                    int l ///<angular momentum  
                   ){
//c *** sigmal(eta,l)= Coulomb phase = arg[Gamma(l+1+i*eta)]
//c *** all Coulomb phase from l=0 to given l are placed in slc(l)
//c *** so that sigma_l=slc(l+1) (in radians)
      double C=0.577215664901532860606512; //can use euler constant.
      double  eps=1.0E-13;
      double sigma=-C*eta;
      double rl=double(l);
      double rm=0.0;
      double c1m,cm;
   l2:  
      c1m=eta/(1.0+rm);
      cm=c1m-atan(c1m);
      if(fabs(cm)<eps) goto l1;
      sigma=sigma+cm;
      rm=rm+1.0;
      goto l2;
   l1: 
      for (int m=1;m<=l;m++)
        sigma+=atan(eta/double(m));
      return sigma;
}
	
/** \brief Compute Coulomb Phase for a given charge, mass and energy
*/
double CoulombPhase(double mu,///< Reduced mass  
						int charge, ///< product of charges 
						double energy, ///< CM energy 
						int l ///< Angular momentum
	) {
		double eta=phy::alpha*charge*sqrt(mu/2.0/energy);
		return mav::CoulombPhase(eta, l);
}
  
double Gamma(int n)
{
  double result=1.0;
  if (n<0) FatalError("Error in Gamma function: negative argument");
  for (int i=1;i<n;i++) result*=double(i);
  return result;
}

double Gamma(double x)
{
  if (x==int(x)) return Gamma(int(x));
  if ((x<1.0)&&(x>0.0)) {
    /*using series for 1/Gamma(x) 2 vol Principa Press, Bloomington Ind,1933*/
    double c[26];
    c[0]=1.0;
    c[1]= 0.5772156649015329;
    c[2]=-0.6558780715202538;
    c[3]=-0.0420026350340952;
    c[4]= 0.1665386113822915;
    c[5]=-0.0421977345555443;
    c[6]=-0.0096219715278770;
    c[7]= 0.0072189432466630;
    c[8]=-0.0011651675918591;
    c[9]=-0.0002152416741149;
    c[10]= 0.0001280502823882;
    c[11]=-0.0000201348547807;
    c[12]=-0.0000012504934821;
    c[13]= 0.0000011330272320;
    c[14]=-0.0000002056338417;
    c[15]= 0.0000000061160950;
    c[16]= 0.0000000050020075;
    c[17]=-0.0000000011812746;
    c[18]= 0.0000000001043427;
    c[19]= 0.0000000000077823;
    c[20]=-0.0000000000036968;
    c[21]= 0.0000000000005100;
    c[22]=-0.0000000000000206;
    c[23]=-0.0000000000000054;
    c[24]= 0.0000000000000014;
    c[25]= 0.0000000000000001;

    double sum=0;
    double y=x;
    for (int k=0;k<26;k++) {
      sum+=c[k]*y;
      y*=x;
    }
    return 1/sum;
  }

  if (x>=1) {
    double y=x-int(x);
    double product=1.0;
    int k;
    for(k=0;k<int(x);k++) product*=y+k;
    return product*Gamma(y);
  }

  return PI/(Gamma(1-x)*sin(PI*(1-x)));
  /* using Gamma(x)*Gamma(1-x)=PI*csc(PI*x) */
}

double Beta(int x,int y)
{
  return Gamma(x)*Gamma(y)/Gamma(x+y);
}

double GarmonicSeries(int n)
{
  double result=0.0;
  for (int i=1;i<n+1;i++) result+=1.0/i;
  return result;
}

double LogGamma
(
    double x    // x must be positive
)
{
	if (x <= 0.0) CriticalError( "Invalid input argument " << x <<  ". Argument must be positive.") ;
	

    if (x < 12.0)
    {
        return log(fabs(Gamma(x)));
    }

	// Abramowitz and Stegun 6.1.41
    // Asymptotic series should be good to at least 11 or 12 figures
    // For error analysis, see Whittiker and Watson
    // A Course in Modern Analysis (1927), page 252

    static const double c[8] =
    {
		 1.0/12.0,
		-1.0/360.0,
		1.0/1260.0,
		-1.0/1680.0,
		1.0/1188.0,
		-691.0/360360.0,
		1.0/156.0,
		-3617.0/122400.0
    };
    double z = 1.0/(x*x);
    double sum = c[7];
    for (int i=6; i >= 0; i--)
    {
        sum *= z;
        sum += c[i];
    }
    double series = sum/x;

    static const double halfLogTwoPi = 0.91893853320467274178032973640562;
    double logGamma = (x - 0.5)*log(x) - x + halfLogTwoPi + series;    
	return logGamma;
}

}


#endif  // __GAMMA__ //

