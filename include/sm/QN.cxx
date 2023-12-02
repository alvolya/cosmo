/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file QN.cxx
\brief Helpful functions for Quantum Numbers
*/
/*
03/14/2004 Separated QN from System
06/29/2009 Made Symmetry operator a part of this file
 */
#ifndef __QN__CXX__
#define __QN__CXX__
#include <tav/basic.cxx>
using tav::Round;
#include <cmath>
#include <string>
#include "spinclass.cxx" //if there is no spinclass just use typedef int spin;
#include "parityclass.cxx" //if there is no spinclass just use typedef int spin;

namespace SM {
  /*!Define consept of spin which is 2*l
\todo make that an equivalent class, typedef still does not allow operator override
*/
  //typedef int spin;
/*! Test if double is an integer or half integer
*/
  bool SpinTest(double x, double prec=0.1) { 
    if (fabs(2.0*x-Round(2.0*x))>prec) return false;
    else return true;
  }
/*! Extract spin from \$ \langle J^-J^+\rangle =J(J+1)-M(M+1) \$
\return Function returns spin if ok or  -1 (unphysical spin -1/2) if the result is unphysical. 
*/
  spin JJSpin(
double x, ///< Input  \$ \langle J^-J^+ \rangle\$
  spin M, ///< Input magnetic projection
double prec=0.1 ///< Input tolerance in deviation from (half)integer
) { 
   double spinout=((sqrt(4*x+double(int(M)*(int(M)+2))+1.0)-1.0)); //twice spin 
//cerr<<double(M)<<" "<<M<<endl;
 //  double spinout=((sqrt(4.0*(x+double(M)*(double(M)+1.))+1.0)-1.0));//double twice spin 
    //cerr<<" spin= "<<spinout<<endl;
    //spin out put is 2*J=spin must be integer
    //now test it
	double testint=0.5*(spinout+double(int(M))); //this is J+M must be integer
	//cerr<<"spintest "<<testint<<endl;
    if (fabs(testint-double(tav::Round(testint)))<prec) 
      return tav::Round(spinout);
    else {
#ifndef COSMONDEBUG
        std::cerr<<"Bad spin: "<<spinout/2.0<<" J-J+= "<<x<<'\n';
#endif
        return -1; //return fail}
    }
  }

  //! set of abelian quantum numbers, additive
  struct AbelianQN {
    //!parity
    int P; //parity
    //!projection Jz
    spin Jz; 
    //!isospin projection
    spin Tz; //projection Tz
#ifdef LScoupling
      int Lz;
      SM::spin Sz;
#endif
  };

    
    AbelianQN & operator+=(AbelianQN &x, const AbelianQN &y) {
        x.P^=y.P;
        x.Jz+=y.Jz;
        x.Tz+=y.Tz;
#ifdef LScoupling
        x.Lz += y.Lz;
        x.Sz += y.Sz;
#endif
        return x;
    }
    
  AbelianQN operator + (const AbelianQN &x, const AbelianQN &y)
    {
        AbelianQN p;
        p.P = x.P^y.P;
        p.Jz = x.Jz + y.Jz;
        p.Tz = x.Tz + y.Tz;
#ifdef LScoupling
        p.Lz = x.Lz + y.Lz;
        p.Sz = x.Sz + y.Sz;
#endif
        return p;
  };

    AbelianQN operator - (const AbelianQN &x, const AbelianQN &y) {
    AbelianQN p;
    p.P=x.P^y.P;
    p.Jz=x.Jz-y.Jz;
    p.Tz=x.Tz-y.Tz;
#ifdef LScoupling
        p.Lz = x.Lz - y.Lz;
        p.Sz = x.Sz - y.Sz;
#endif
    return p;
  };
    
    bool operator == (const AbelianQN &x, const AbelianQN &y) {
#ifdef LScoupling
        return ((!(x.P^y.P))&&(x.Jz==y.Jz)&&(x.Tz==y.Tz)&&(x.Lz==y.Lz)&&(x.Sz==y.Sz))
#endif
        return ((!(x.P^y.P))&&(x.Jz==y.Jz)&&(x.Tz==y.Tz));
    }
    
    bool operator != (const AbelianQN &x, const AbelianQN &y) {
        return (!(x==y));
    }
    
   //! set of abelian quantum numbers, additive
  struct NonAbelianQN {
         //!particle number
    int N; //particle number
    //! total spin
    spin J;
    //! isospin
    spin T;
    int L; ///< orbital angular momentum
	spin S; ///<spin quantum number 1/2
	 std::string tag; ///<label
  };
  

  
  bool operator == (const NonAbelianQN &a, const NonAbelianQN &b) {
  return ((a.N==b.N)&&(a.J==b.J)&&(a.T==b.T)&&(a.L==b.L)&&(a.S==b.S));
  }
  
  bool operator != (const NonAbelianQN &a, const NonAbelianQN &b) {
  return (!(a==b));
  }

  //!\brief Set of quantum numbers
  struct QN: AbelianQN, NonAbelianQN{
  };


 bool operator == (const QN &a, const QN &b) {
  return ((AbelianQN(a)==AbelianQN(b))&&(NonAbelianQN(a)==NonAbelianQN(b)));
  }
   bool operator != (const QN &a, const QN &b) {
  return (!(a==b));
  }

  

///Print Quantum numbers  
std::ostream& operator << (std::ostream& stream, const QN q)
{
//  stream<<q.N<<'\t';
//  stream<<q.J<<'\t';
//  stream<<q.Jz<<'\t';
//  stream<<q.T<<'\t';
//  stream<<q.Tz<<'\t';
//  stream<<q.L<<'\t';
//  stream<<q.P<<'\t';
    stream<<"N="<<q.N<<' ';
    stream<<"J="<<q.J<<' ';
    stream<<"Jx="<<q.Jz<<' ';
    stream<<"T="<<q.T<<' ';
    stream<<"Tz="<<q.Tz<<' ';
    stream<<"L="<<q.L<<' ';
    stream<<"P="<<q.P<<' ';
        return stream;
}

///Read quantum numbers 
template <class cnumber>  
std::istream& operator >> (std::istream& stream, QN &q)
{
  stream>>q.N;
  stream>>q.J;
  stream>>q.Jz;
  stream>>q.T;
  stream>>q.Tz;
  stream>>q.L;
  stream>>q.P;
        return stream;
}
/*!\brief Operator with symmetry. Tags quantum numbers (SM::QN) to any structure
  */
template <typename coperator>
class SymmetryOperator : public coperator, public QN {
  public:
  SymmetryOperator () {}; ///< Default Constructor with nothing.
  SymmetryOperator (const SymmetryOperator &SX):coperator(SX),QN(SX) {}; ///< Constructor from itself. 
  SymmetryOperator (const QN& J, const coperator &T):coperator(T),QN(J) {}; ///<Constructor from pair
  /*! \brief basic constructor from many-body coperator
  \note This will be called if constructor SymmetryOperator(int) is called
  */
  SymmetryOperator (const coperator &T):coperator(T) {}; 
  ///\brief Copy from quantum numbers
  SymmetryOperator& operator = (const QN J){QN::operator=(J); return *this;};
  ///\brief Copy from operator
  SymmetryOperator& operator = (const coperator J){coperator::operator=(J); return *this;};
  //using coperator::operator = ;
 // using QN::operator =;
    /*! Print function for symmetry operator
     */
    
    friend std::ostream& operator<<(std::ostream& os, SymmetryOperator& a) {
        os<<QN(a)<< std::endl;
        os<<coperator(a);
        return os;
    }
};



}
#endif 

