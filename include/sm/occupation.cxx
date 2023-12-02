/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file occupation.cxx
 \brief Computational of diagonal occupation numbers
 
 Using operators to comute occupation number is ineffective. No phase treatment is needed
 
\ingroup gp_sm
*/
#ifndef __OCCUPATION__CXX__
#define __OCCUPATION__CXX__

#include <debug.h>
#include <sm/ManyBodyStates.cxx>
//#include "StatusBar.cxx"
namespace SM {
/** \brief Compute diagonal density for any occupation mapping 
 basically we count how often a certain partition is encountered where partition is defined by the map
*/
template <typename coperator>
int occupation(
coperator &occ, ///<output s.p. density  
double *af, ///< state   
const Many_Body_States &st, ///< reference basis states
int *map //< mapping of single-particle state to level, e.g. map[s.p state]=level 
)
  //using vectors af and ai we compute overlap into matrix ovl[1][2] 
  {
    //we start by making our transition matrix zero
	occ.clear();
	int mapdim=occ.rank; //number of levels in map
    typedef typename coperator::value_type::first_type myinteger;
	myinteger element = new  int [mapdim]; //make temporary storage element
 for (int s1=0;s1<st.n;s1++) 
       {
       for (int ll=0;ll<mapdim;ll++) element[ll]=0; //make element zero
       for (int ni=0;ni<st.N;ni++) element[map[st.z[s1][ni]]]++; //find our occupation distribution for this mapping       
	//incert this distribution into map
	   occ[element]+=tav::Sqr(std::abs(af[s1]));	  
     } //end state scanning
delete [] element;
      return 0;
  } //end function

    /** \brief This is a very basic occupation numbers
     */
    template <typename coperator>
    int occupation(
                   coperator &occ, ///<output s.p. density
                   double *af, ///< state
                   const Many_Body_States &st ///< reference basis states
                   )
    //using vectors af and ai we compute overlap into matrix ovl[1][2]
    {
        //we start by making our transition matrix zero
        for (spsint i=0;i<occ.size();i++) occ[i]=0;
        for (uint_nbasis s1=0;s1<st.n;s1++)
        {
            for (int ni=0;ni<st.N;ni++) occ[st.z[s1][ni]]+=tav::Sqr(std::abs(af[s1]));	  
        } //end state scanning
        return 0;
    } //end function

  
} //end of name space SM


#endif //__OCCUPATION__CXX__

