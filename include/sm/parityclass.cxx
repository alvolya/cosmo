/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


#ifndef __PARITYCLASS__CXX__
#define __PARITYCLASS__CXX__
#include <debug.h>
//#include <string>
namespace SM{
/*!
\class parity 
\brief parity variable equivalent to boolean
integer conversion is defined as mathematical parity
*/
class  parity {
  private:
  bool bp;
  public: 
    ///default constructor, does nothing
    parity(void) {bp=false;};
    parity(const int i){bp=i&1;};
	parity(const char mychar) {
    if (mychar=='-') {bp=true; return;}
	if (mychar=='+') {bp=false; return;}
    CriticalError("Use '+' pr '-' for parity, symbol '"<<mychar<<"' is unclear ");
	}
	parity(const char *mystring) {
	*this=mystring[0];
    //bp=(mystring[0]=='-') 
	}  
	//parity(const std::string mystring) {*this=mystring.c_str(); }
	

 /// Elegant output 
    friend std::ostream& operator << (std::ostream& stream, const parity &s) {
    stream<<char(s); return stream;};
    friend std::istream& operator >> (std::istream &stream, parity &s) {
	char cc; stream>>cc; s=cc; return stream;
    };


    //This is the most important operator making parity equivalent to integer in implicit conversions
	/*! \brief Set equivalence to boolean -/+ for true/false */
	inline operator bool(void) const {return bp;};
	/*! \brief Set equivalence to character + or -
	*/
	operator char(void) const {if (bp)  return '-'; else return '+';}
	
	/*! \brief Parity multiplication, composite system 
	*/
	parity operator*=(const parity a) {bp^=a.bp; return *this;};
	parity operator*(const parity a) {return (this->bp)^a.bp;};
	
	/*! \brief Apply parity (multiply by -1 or +1) with operator () 
	*/
	template <typename cnumber> 
	cnumber operator() (const cnumber x) {return ((bp)? -x :x); };	
 } ; //end class
 

 
} //end SM namespace

#endif //end __PARITYCLASS__CXX__

