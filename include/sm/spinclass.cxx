/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*!
\author Alexander Volya  <http://www.volya.net>
\class spin
\brief half-integer spin varable
Class includes implicit conversion to double and construction from
integer as well as double;
- Variable fully behaves as integer 2i=s
- Elegant read
- include string class for easier print
Note we need to choose a phylosophy, Chosen (a)
spin is always implicitly converded to integer and we deal with that
We do not define + because if both spins we still convert to integers, faster internal processing
Rules
 spin=double //will work as expected
 spin=int       //will assume spin eqals to 2*int
 result of any mathematical operaton implicity converst to integer and result is an integer
 integers must be implicitly cast to spin. 

*/
#ifndef __SPINCLASS__CXX__
#define __SPINCLASS__CXX__
#include <cstdio>
using std::sprintf;
#include <cstdlib>
#include <cstring>
#include <cmath> //for std::fabs
#include <iostream>
using std::atoi;
using std::atof;
using std::strcpy;
using std::strchr;
#include <string>
//using namespace std;
#define FatalError_spinclass(message)              { std::cerr <<"Error: "<< message <<std::endl<< __FILE__ << ":" << __LINE__ << std::endl; throw (1); }
namespace SM{

class  spin {
  private:
  int twicespin;
  public: 
    ///default constructor, does nothing
    spin(void) {};
    ///constructor from integer
    spin(const int i):twicespin(i) {};
    /// constructor from itself 
    spin(const spin &a):twicespin(a.twicespin) {};
    /// constructor from double
    spin(double x) {
        if (x>=0) twicespin=(int)(2.0*x+0.5);
        else twicespin=(int)(2.0*x-0.5);
      if (std::fabs(twicespin-2.0*x)>0.01) FatalError_spinclass("spin: constructor spin(double): spin "<<x<<" is not integer or half-integer \n");
    } 
	///constructor from cstring
	spin(const char *mystring) {
	const int ssize=16;//maximum input size
	char ctmp[ssize]; //temporary char
    std::strcpy(ctmp,mystring);
	char *pch=std::strchr(ctmp,'/'); //develop later pch==NULL not found
	  if (pch!=NULL) { //found slash
	  if (pch[1]!='2') FatalError_spinclass("Error: operator spin(const char*): bad spin "<<ctmp);
	  (pch[0])='\0'; //replace / with termination
	  twicespin=std::atoi(ctmp);
	  }
	 else {
//this is double: note it maybe junk but we can not control everything
	*this=std::atof(ctmp); //here we could use an already avalable function
	//twicespin=Round(2.0*x);
     // if (tav::Abs(twicespin-2.0*x)>0.01) FatalError("spin: constructor spin(double): spin "<<x<<" is not integer or half-integer \n");
	}

    } 
	
		
	spin(const std::string mystring) {*this=mystring.c_str(); }
	
    char * SpinToChar(char *mystring) {
		const int ssize=16;//maximum input size
    if (twicespin&1)  snprintf (mystring,ssize, "%i/2", twicespin); 
      else snprintf (mystring,ssize, "%i", (twicespin>>1)); 
	  return mystring;
	}
    /// Elegant output 
    friend std::ostream& operator << (std::ostream& stream, const spin &s) {
      if (s.twicespin&1) stream<<s.twicespin<<"/2";
      else stream<<(s.twicespin>>1);
      return stream;
    };
    friend std::istream& operator >> (std::istream &stream, spin &s) {
      const int ssize=16;//maximum input size
      char ctmp[ssize]; //temporary char
//      const char terminate='\0'; //define terminating character
      stream>>ctmp;
	  s=ctmp; //use string char conversion
	/*  
	  char *pch=std::strchr(ctmp,'/'); //develop later pch==NULL not found
	  if (pch!=NULL) { //found slash
	  if (pch[1]!='2') FatalError("Error: operator << spin: bad spin "<<ctmp);
	  (pch[0])='\0'; //replace / with termination
	  s=std::atoi(ctmp);
	  }
	 else {
//this is double: note it maybe junk but we can not control everything
	s=std::atof(ctmp);
	}
	*/
	
      return stream;
    };

    //This is the most important operator making spin equivalent to integer in implicit conversions
	inline operator int(void) const {return twicespin;};
	//this operator returns a string
	operator std::string(void) const {
	const int ssize=16;//maximum input size
	char mystring[ssize];
	if (twicespin&1)  snprintf (mystring,ssize, "%i/2", twicespin); 
      else snprintf (mystring,ssize, "%i", (twicespin>>1)); 
	  return std::string(mystring);
	}
	
	spin operator+=(const spin a) {twicespin+=a.twicespin; return *this;};
	spin operator-=(const spin a) {twicespin-=a.twicespin; return *this;};
		//int operator=(const spin a); 
 } ; //end class
 

 
} //end SM namespace

#endif //end __SPINCLASS__CXX__

