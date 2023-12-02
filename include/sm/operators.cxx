/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file operators.cxx
\brief s.p. operators
\date 10/19/2007


\par Changes
- 11/26/06 Phase in Annihilation operator is based on nnf.
- 07/22/07 Annihilation of more particles than in state is fixed, may not be
         very efficient, check later
- 10/16/07 Operator1B one-body operator

\example operators_test.cpp
\example operators_Creation_Annihilation_test.cpp
*/

/*! \page CI M-scheme operators and states
\section sec_phase Phase Conventions
-#  Boolean variable phase is true, when phase is needed:
          -# True:  phase=-1 odd C++ True=1
          -# False: phase=+1 even C++ False=0
-#  Summary of boolean operations and related parity performed by XOR "^"
 0^0=0, 1^1=0 but 1^0=1 and 0^1=1
        -#    -1*(-1)=1 odd*odd=even       1^1=0     
        -#              odd*even=odd     1^0=1     
        -#     -1*1=-1  even*odd=odd     0^1=1     
        -#     1*1=1   even*even=even    0^0=0
-#  <var>n&1</var> is phase of \f$(-1)^n\f$: i.e.  1 (True) for odd number, 0(False) for even 
-#  to apply phase  <var>result= (phase? -x : x )</var>;
\section sec_operator m-scheme many-body state/operator
An many-body state/operator in the m-scheme 
is coded by an array SM::spsint *a[]={0,1,2..} 
which is interpreted as an operator \f$a^\dagger_0 a^\dagger_1 \dots \f$. 
Ordering is assumed 0<1<2... . The m-scheme 
ket-state is \f$|\alpha \rangle = a^\dagger_0 a^\dagger_1 \dots|0\rangle \f$. 
The number of s.p. operators in a state or many-body operator must be specified separately. 

\section sec_cioperator many-body operator 
 The full many-body operator is reprented by the Sparse Tensor tav::STensor which is 
 a sum of \ref sec_operator with corresponding coefficients.
 An example of a triplet creation operator can be \f$ \hat{T}^\dagger=\sum_{123}\, C_{123} a^\dagger_1 a^\dagger_2 a^\dagger_3 \f$ which is coded as a 
 set {123} mapped to \f$ C_{123}\f$, call it <var> rank={0+3) </var> . However, the fact that this is a creation operator must be identified and  all possible interpretation of the structure are 
- \f$ \hat{T}^\dagger=\sum_{123}\, C_{123} a^\dagger_1 a^\dagger_2 a^\dagger_3 \f$ rank=(0+3)
- \f$ \hat{T}=\sum_{123}\, C_{123} a_3 a_2 a_1 \f$ rank=(3+0)
- \f$ \hat{X}=\sum_{123}\, C_{123} a^\dagger_2 a^\dagger_3 a_1 \f$ rank=(1+2)
- \f$ \hat{Y}=\sum_{123}\, C_{123} a^\dagger_3 a_2 a_1 \f$ rank=(2+1).

 Thus, a sparse tensor mapping must be suplimented with explicit number of annihilation (first part of the set) and creation operators (second part of the set) in it. The normal ordering form is always assumed. 
\section sec_forward forward ordering of Hermitian operators 
Hermitian operators \f$ H^\dagger=H \f$ in m-scheme basis are stored in upper-right matrix form, which we refer to as forward ordering. The forward ordered operator 
\f$\acute{H}|i\rangle=|f\rangle \f$ where in m-scheme ordering \f$ i\le f \f$.
*/ 

#ifndef OPERATORS__CXX__
#define OPERATORS__CXX__
//#include <tav/sort.cxx> 
//using tav::Locate;
#include <cstring>
using std::memcpy;	
/*! \namespace SM 
\brief Shell Model / Many-body Configuration Interaction  
*/

namespace SM{
  
/*! \brief Many-body comparison 
Determine the minimum number of s.p. operations needed to make final state from initial
 */
  int Compare(
  int Nopr, ///< number of s.p. operators in operator
  spsint *bra, ///< operator 1
  spsint *ket ///< operator 2
                   )
  {
    int rank=Nopr;
    int i=0;
    int j=0;
    for (i=0;i<Nopr;i++) {
    for (;(j<Nopr)&&(bra[i]>ket[j]);j++) ; 
    if (bra[i]==ket[j]) rank--;
    }
    
  return rank;
  }
  
/*! \brief Many-body comparison 
Determine the minimum number of s.p. operations needed to make final state from initial the states are in arbitrary order
 */
  int CompareGeneral(
  int Nopr, ///< number of s.p. operators in operator
  spsint *bra, ///< operator 1
  spsint *ket ///< operator 2
                   )
  {
    int rank=Nopr;
    for (int i=0;i<Nopr;i++) for (int j=0;j<Nopr;j++) if (bra[i]==ket[j]) rank--;
  return rank;
  }  

/*! \brief Many-body Annihilation operator 

This is annihilation operator in a string format \sa \ref sec_operator
- The state {1,2,3,4...} \f$ \rightarrow a^\dagger_1 a^\dagger_2 a^\dagger_3 a^\dagger_4
\dots |0\rangle \f$ 
- operator {1,3} \f$ \rightarrow a^\dagger_1 a^\dagger_3 \f$ and
as annihilation is \f$ a_3 a_1 \f$ 
- result \f$ -a^\dagger_2 a^\dagger_4 \dots |0\rangle \f$ interpreted as {2,4,...} and phase=true
\note this operator is faster then version SM::Annihilation1() 
\return The boolean variable is returned False if action can not be done


 */
    bool Annihilation(
bool &phase, ///< phase change is applied to this variable, \sa \ref sec_phase
int &Nket, ///< number of s.p. operators in final state 
spsint *ket, ///< the final state 
int Nopr, ///< number of s.p. operators in operator
spsint *opr, ///< operator
int Nbra, ///< number of s.p. operators in initial state
spsint *bra ///< initial state vector
)
{
  if (Nopr>Nbra) return false; //fix problem with one particle, review later
      int nnf=0;
      int nni=0;
      int nno=0;
      int itmp;
      
      for (;nno<Nopr;nno++) //scan through operator 
	{  itmp=(Nbra-nni)-(Nopr-nno);
	  for (;bra[nni]<opr[nno];nni++) 
	    if (itmp--) {ket[nnf]=bra[nni]; nnf++; } 
	    else return false;  //just copy unless we need to exit
	  //Nbra-nni=Nopr-nno
	  if  (opr[nno]!=bra[nni]) return false; 
	  else {nni++; phase^=(nnf&1);} // move forward  
      }
      //need to finish up
      //std::cout<<"finishup \n";
      if (nni-nnf-Nopr)  return false;
      else {
	for (;nni<Nbra;nni++) {ket[nnf]=bra[nni]; nnf++;}}
      Nket=Nbra-Nopr; //number of particles in a final state, compute at end
      return true;
}


/*! \brief Maby-body Creation operator

See notes SM::Annihilation()

 */
bool Creation(
    bool &phase, ///<phase to be applied
    int &Nket, ///<size of final vector
    spsint *ket, ///<final vector
    int Nopr, ///<size of operator
    spsint *opr, ///<operator
    int Nbra, ///< size of initial state 
    spsint *bra ///<size of final state 
)
{
      int nnf=0;
      int nni=0;
      int nno=0;
      Nket=Nbra+Nopr; //number of particles in a final state
      while (nnf<Nket) //scan through state
	{ 
	  if (nni==Nbra) {ket[nnf]=opr[nno];  phase^=(nni&1);// add phase 
	    nnf++; nno++; continue;}  //can continue with operator only
	  if (nno==Nopr) {ket[nnf]=bra[nni]; nnf++; nni++; continue;} 
	  if (bra[nni]<opr[nno])  {ket[nnf]=bra[nni]; nnf++; nni++; continue;}
	  if (bra[nni]>opr[nno]) {ket[nnf]=opr[nno];  phase^=(nni&1);// add phase 
	    nnf++; nno++; continue;}
	  return false;
	}

      return true; 
}

/// @cond DEV
/*
This is annihilation operator opr|bra>=|ket> phase
everything is ordered bra[i]<bra[i+1]
 */
bool Annihilation1(bool &phase, 
int &Nket, spsint *ket, int Nopr, spsint *opr, int Nbra, spsint *bra )
{
      int nnf=0;
      int nni=0;
      int nno=0;
      Nket=Nbra-Nopr; //number of particles in a final state
      //fully copy or annihilate initial state
      for (;nni<Nbra;nni++) //scan through state
	{ 
	  
	  if ((bra[nni]==opr[nno])&&(nno<Nopr)) {//a particle can be removed
	    //   phase^=(nni&1); //phase in ket
	    // phase^=(nno&1); //phase in operator
	    phase^=(nnf&1); //phase on difference
	    nno++; //this will move to the next particle in operator
	    
	  }
	  else if (nnf<Nket) {ket[nnf]=bra[nni]; nnf++;} else return false;
	  //nni++; //always shift initial state

	}
      // if (nno!=Nopr) return false; //this is a falure answer is zero;
      //else return true;
      return true;
}

bool Creation1(bool &phase, 
int &Nket, spsint *ket, int Nopr, spsint *opr, int Nbra, spsint *bra )
{
      int nnf=0;
      int nni=0;
      int nno=0;
      Nket=Nbra+Nopr; //number of particles in a final state
      for (;nno<Nopr;nno++) {
	for (;(bra[nni]<opr[nno]);nni++)  {
	  if (nni==Nbra) { //just finish 
	    for (;nno<Nopr;nno++) {ket[nnf]=opr[nno]; nnf++; phase^=(nni&1); 
	     }  return true;
	  }
	  ket[nnf]=bra[nni]; nnf++; } 
        //exit loop 
	if (bra[nni]==opr[nno]) return false; //fermions exit
	else {ket[nnf]=opr[nno]; nnf++; phase^=(nni&1); } //put particle
      }
      //operator is out just copy
      for (;nni<Nbra;nni++)  {
	  ket[nnf]=bra[nni]; nnf++; } 
            return true; 
}
/// @endcond

//! logarithm base 2 of SM::spsint used as x*sizeof(spsint)=x<<log_char_size
//#define log_char_size 0
//this is an extremely dangerous place log_char_size is log_2(sizeof(spsint))
//2=int 0=[unsigned] char

// /*! \brief One-body operator 
// This is operator of the form \f$a^\dagger_2 a_1 \f$ given as array {12} 
// the state vector is ordered {1234..N} and corresponds to \f$a^\dagger_1 a^\dagger_2 \dots a_N^\dagger\f$ where  1<2<3<... 
// \return The boolean variable is returned False if action can not be done
// \note It is faster to copy an entire array compared to copy by parts.
//  */
//     bool Operator1B(
//     bool &phase, ///< phase change is applied to this variable, \sa \ref sec_phase
//     spsint *ket, ///< final state vector
//     const spsint *opr, ///< operator opr[1] for opr[1]=f(to be created) opr[0]=i (annihilated) 
//     int NN, ///< number of s.p. operators in initial state
//     const spsint *bra ///< initial state vector
// 		     )
// {	
// 	int i=tav::Locate(NN,bra,opr[0]); //locate particle to be removed
// 	if (i<0) return false; //operation can not be performed, no particle
// 	if (bra[i]!=opr[0]) return false; //opration can not be performed
// 	if (opr[1]==opr[0]) {
// 		memcpy(ket,bra, NN<<log_char_size);
// 		//for (int ii=0;ii<NN;ii++) ket[ii]=bra[ii]; //copy an array
// 		return true; 
// 	}
// 	int f=tav::Locate(NN,bra,opr[1]); //locate particle to be created
// 	//order is ascending so
// 	if (bra[f]==opr[1]) return false; //opration can not be performed
// 	if (i<=f) {
// 		//we can alternatively copy an entire array
// 		memcpy(ket,bra,NN<<log_char_size); 
// 		//** memcpy(ket,bra,i<<log_char_size); //copy all elements [0,i)
// 		memcpy(ket+i,bra+i+1,(f-i)<<log_char_size); 
// 		//copy ket[i,f)=bra[i+1,f] 
// 		ket[f]=opr[1]; //incert ket
// 		//** memcpy(ket+f+1,bra+f+1, (NN-f-1)<<log_char_size);
		
// 		phase^=((i+f)&1);
// 	} 
// 	else {  //below line commented -> uncomment  all //**
// 		memcpy(ket,bra,NN<<log_char_size); 
// 		f++; //make f->f+1
// 		//** memcpy(ket,bra,f<<log_char_size); //ket [0,f+1)= bra[0,f+1)
// 		//----------core operations--------------
// 		ket[f]=opr[1]; //ket[f+1]=operator
// 		memcpy(ket+f+1,bra+f,(i-f)<<log_char_size); 
// 		//ket [f+2,i+1) =bra[f+1,i)
// 		//---------end core operation
// 		phase^=((i+f)&1);
// 		//**  i++; //make i->i+1 to speed what is below
// 		//**  memcpy(ket+i,bra+i, (NN-i)<<log_char_size);
// 		 //ket[i+1,N)=bra[i+1,N)
		
// 	}
// 		return true;
// }

}
#endif 


