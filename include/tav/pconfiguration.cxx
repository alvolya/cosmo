/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file pconfiguration.cxx
\brief Particle configuraitons
\ingroup gp_partitions
  
Configuration in the pp format is represented by an array of locations
The locations are ordered in increasing position. 
The order of configurations is "inverted"
 
*/
#ifndef __PCONFIGURATION__CXX__
#define __PCONFIGURATION__CXX__

namespace tav {

double  fBinomialCoefficient(int W, int N) {
if (N>W) return 0.0;
double C = 1.0;
for (int kk = W; kk > W-N; kk--) C*=double(kk);
for (int kk=1;kk<=N;kk++) C/=double(kk);
return C;
}


//double version computing fraction
template <typename vspsint>
double ConfigurationPosition(
int W, ///< number of single particle states
int  N, ///< Number of particles
vspsint &currentstate ///< Current operator - input; Next state - output
) {
double q=fBinomialCoefficient(W,N);
double f=q;
for (int i=0;i<N;i++) {q-=fBinomialCoefficient(W-currentstate[i]-1,N-i);
}
return q;
}

//double version computing fraction
template <typename vspsint>
double FractionalConfigurationPosition(
int W, ///< number of single particle states
int  N, ///< Number of particles
vspsint &currentstate ///< Current operator - input; Next state - output
) {
double q=fBinomialCoefficient(W,N);
double f=q;
for (int i=0;i<N;i++) {q-=fBinomialCoefficient(W-currentstate[i]-1,N-i);
}
return q/f;
}

/*! Generate Next Fermi Distribution 
Given a many-body state (see. \ref CI ) we generate a next Fermi-allowed 
distribution. Next corresponds to moving last particle forward
particles = 3  s.p. space dimension 5


dist    repr         state

--------------------------------

xxx..   {012}        state 0

xx.x.   {013}        state 1

xx..x   {014}         state 2

x.xx.   {023}         state 3

\return Boolean flag true if successful, false if not (reached end).
*/
template <typename vspsint>
bool NextFermiDistribution(
int W, ///< number of single particle states
int  ph, ///< Number of particles 
vspsint &currentstate ///< Current operator - input; Next state - output
) {
  //currentstate has dimention ph 
  //there are ph things and they can take values from 0 to N-1
      //iterations are here
  if (ph==0) return false; //if ph==0 return false 
    bool  OKflag=false;
    for (int i=ph-1;i>=0;i--) {
  //go in reversed order
      /* I have ph-i particles moving starting from I 
      [x-xx-x-Xxx-x-x]
       0 1    i     W-1                 
      */    
      //  if (currentstate[i]+1==N) continue; //we are at end [x-xxx----X]
      // if ((i!=(ph-1))&&((currentstate[i]+1)==currentstate[i+1])) continue;
      // we need to move i-th particle and dump all others after it
      // [x--xX---xxx]->[x--x-Xxxx---]
      //number of spaces I have after it is W-1-currentstate[i]
      //number of particle including X is ph-i
      if ((W-1-currentstate[i])<(ph-i)) continue; //we do not fit;
      //there is a trick here we need to change current i last
      for (int j=ph-1;j>=i;j--) {currentstate[j]=currentstate[i]+j+1-i; }
      // n++;


      OKflag=true;
      break;
    }
    return OKflag;
}
/*! Generate Previous Fermi distribution, see NextFermiDistribution()
*/
template <typename vspsint>
bool PreviousFermiDistribution(
int W, ///< number of single particle states
int  ph, ///< Number of particles
vspsint &currentstate ///< Current operator - input; Previous state - output
) {
  //currentstate has dimention ph 
  //there are ph things and they can take values from 0 to W-1
      //iterations are here
 if (ph==0) return false; //if ph==0 return false 
     int i=ph-1;
    while ((i>0)&&(currentstate[i-1]+1==currentstate[i])) i--; 
    //at the end i will tell end seriens [x--xxxx---]
    //                                       i  ph-1
    //now the i tells me how many I have in a row if i=ph-1 then non
    //if i=0 then all are in a row
    if ((i==0)&&(currentstate[i]==0)) return false;
    //now we are sure that we have no
    currentstate[i]--; //move one step back last particle
    int in=W-1; //set to the end
    // for (int j=i+1;j<ph;j++) {currentstate[j]=in; in--;}
    for (int j=ph-1;j>i;j--) {currentstate[j]=in; in--;}
  //go in reversed order same as with Next
     
    return true;
}

/*! Generate Next Boson distribution, see NextFermiDistribution()
*/
template <typename vspsint>
bool NextBosonDistribution(
int W, ///< number of single particle states
int  ph, ///< Number of particles
vspsint &currentstate ///< Current operator - input; Previous state - output
) {
  //currentstate has dimention ph 
      //iterations are here
 if (ph==0) return false; //if ph==0 return false 
    bool  OKflag=false;
    for (int i=ph-1;i>=0;i--) {
  //go in reversed order
      /* I have ph-i particles moving starting from I 
      [x-xx-x-Xxx-x-x]
       0 1    i     W-1                 
      */    
      //  if (currentstate[i]+1==W) continue; //we are at end [x-xxx----X]
      // if ((i!=(ph-1))&&((currentstate[i]+1)==currentstate[i+1])) continue;
      // we need to move i-th particle and dump all others after it
      // [x--xX---xxx]->[x--x-Xxxx---]
      //number of spaces I have after it is W-1-currentstate[i]
      //number of particle including X is ph-i
      //bosons
      if (currentstate[i]==(W-1)) continue; //we are at end;
      //there is a trick here we need to change current i last
      // for (int j=ph-1;j>=i;j--) {currentstate[j]=currentstate[i]+j+1-i; }
      // n++;
 for (int j=ph-1;j>=i;j--) 	currentstate[j]=currentstate[i]+1;
	//put quantum numbers back

      OKflag=true;
      break;
    }
    return OKflag;
}



}//end namespace tav
#endif //__PARTITION__CXX__

