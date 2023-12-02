/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*!\defgroup gp_mbs Many-body states
Many-body states class and generation of many-body states 
This is a dedicated class to many-body states which links integer number to a particular m-scheme operator that creates that state, see \ref CI .
*/
/*!
\author Alexander Volya  <http://www.volya.net>
\file ManyBodyStates.cxx
\brief Many-body states class
\ingroup gp_mbs
*/
/*
10/29/09 MakeMBS2 bug fix
04/13/2006 Multiple definitions are excluded
03/18/2006 typdef unsigned char spsint is in this file
 */
#ifndef __MANYBODYSTATES__CXX__
#define __MANYBODYSTATES__CXX__
#include <debug.h>
//typedef unsigned char spsint;
//typedef unsigned long long  uint_nbasis;
#include <sm/SingleParticleStates.cxx>
#include <iostream>
using std::ofstream;
using std::ifstream;
using std::ios;
using std::cout;

//#include <sm/System.cxx>
namespace SM {

int compare(int N, spsint *a, spsint *b) ;








/*!\brief Many-body states class
\ingroup gp_mbs
*/
class Many_Body_States {

public:
  int N; ///<number of particles 
  uint_nbasis n; ///<number of states
  spsint **z; ///<states,   z[i]-operator, (see \ref CI ) e.g z[i][n] gives location of n-th particle in i-th state 
  //  Many_Body_States Many_Body_States::operator = (const Many_Body_States &);
  Many_Body_States(void);
  ~Many_Body_States(void);
  int Write(const char *);
  int Read(const char *);
};
///\brief Default constructor 
Many_Body_States::Many_Body_States (void) {
  z=NULL;
  n=0;
}
Many_Body_States::~Many_Body_States (void) {
 if (z!=NULL) {
  for (uint_nbasis nn=0;nn<n;nn++) delete [] z[nn];
  delete [] z;
  z=NULL;
 }
}

/*!\brief compare two operators a?b returns 1 if a>b, 0 for a=b and -1 for a<b
\todo Relate this to ArrayComparison in stensor.cxx
*/
int compare(
int N, ///< number of particles 
spsint *a, ///< operator a
spsint *b  ///< operator b
) {//this will compare two states
  for(int i=0;i<N;i++) {if (a[i]<b[i]) return -1; if (a[i]>b[i]) return 1;}
	//if we get through this and nothing they are equal
	return 0;
  }


    /*! \brief Locate number from a given operator
     \return State number, if locate did not find a state number n+1 (out of range ) is returned.
     */
    uint_nbasis locate(
                       const Many_Body_States &st, ///<Many-body states
                       spsint *a, ///<State-operator to find
                       uint_nbasis jl=0 ///<Search up from lower limit l1
    )
    {
        uint_nbasis jm,ju;
        int ascnd;
        
        //jl=0; //lower limit search
        ju=st.n; //upper limit
        //	cout<<"lower limit "<<jl<<" upper limit "<<ju<<endl;
        int happy=compare(st.N, a, st.z[jl]) ; //happy=1 a>jl
        jm=jl;
        if ((happy*compare(st.N, a, st.z[ju-1]))>0)
        {
            //cerr<<"error in SM::locate :no state in the interval"<<endl;
            return (st.n+1);
        }
        while (happy!=0)
        {
            jm=(ju+jl) >> 1; //computes the mid point (BINARY SHIFT)
            if (jm==jl)
            {
                //cerr<<"state is not in the list "<<endl;
                return (st.n+1);
            }
            //middle point less or = than middle
            happy=compare(st.N, a, st.z[jm]);
            if (happy==1) //a>z[jm]
                jl=jm;
            else
                ju=jm;
            //	        cout<<"lower limit "<<jl<<" upper limit "<<ju<<" mid "<<jm<<" happy="<<happy<<endl;
            
        }
        return jm;
        
    }

/*! \brief Locate number from a given operator, search within range
\return State number, if locate did not find a state number n+1 (out of range ) is returned.
*/
  uint_nbasis locate(
	const Many_Body_States &st, ///<Many-body states 
	spsint *a, ///<State-operator to find
	uint_nbasis jl, ///<lower limit
	uint_nbasis ju ///<upper limit
	) {
        uint_nbasis jm;
        int ascnd;
        ju++; //c notaton increment by one
        //jl=0; //lower limit search
        //ju=st.n; //upper limit
	//	cout<<"lower limit "<<jl<<" upper limit "<<ju<<endl;
	int happy=compare(st.N, a, st.z[jl]) ; //happy=1 a>jl
	jm=jl;
   	if ((happy*compare(st.N, a, st.z[ju-1]))>0) 
      {
	//cerr<<"error in SM::locate :no state in the interval"<<endl; 
	return (st.n+1);}
        while (happy!=0) {
	  jm=(ju+jl) >> 1; //computes the mid point (BINARY SHIFT)
	  if (jm==jl) {
	    //cerr<<"state is not in the list "<<endl; 
	    return (st.n+1);}
	  //middle point less or = than middle
	  happy=compare(st.N, a, st.z[jm]); 
	    if (happy==1) //a>z[jm]
                        jl=jm;
                else
                        ju=jm;
	    //	        cout<<"lower limit "<<jl<<" upper limit "<<ju<<" mid "<<jm<<" happy="<<happy<<endl;
	    
        }
	return jm;

  }


  //locate is much faster......
/*! \brief Search with Hunt method, slower than locate not used
\return State number
*/
  uint_nbasis hunt(const Many_Body_States &st, spsint *a, uint_nbasis jl) 
  {
    uint_nbasis ju=st.n-1;
    uint_nbasis stt=1;
    uint_nbasis previousj=jl;
    uint_nbasis currentj=jl+stt;
    int flag;
    if (locate(st,a,previousj,ju)==0) return previousj;
    while (true) {
    if (currentj>ju) return locate(st,a,previousj,ju); //perform bisection
    flag=compare(st.N, a, st.z[currentj]);
    if (flag==0) return currentj;
    if (flag<0) return locate(st,a,previousj,currentj); //done with search
    previousj=currentj; stt*=2; currentj=jl+stt;
    }
    
  } 




//FILE OPERATIONS
///Save Many-body states to file, write N, n, z
int Many_Body_States::Write (const char *filename) {
  ofstream outf;
  outf.open(filename, ios::binary);
  outf.write((char*)(&N),sizeof(int));
  outf.write((char*)(&n),sizeof(uint_nbasis));
  for (uint_nbasis nn=0;nn<n;nn++) outf.write((char*)(z[nn]), N*sizeof(spsint));
  outf.close();
  return 1;
}
///Read Many-body states from file
int Many_Body_States::Read (const char *filename) {
  ifstream outf;
  outf.open(filename, ios::binary);
if (!outf)   {               // if the file does not exist
  cerr<<"File "<<filename<<" is not found"<<endl; // return error message
  return 0;
}
  //delete old stuff
    if (z!=NULL) {
    cerr<<"Remaking many-body states" <<endl;
    //delete first
  for (uint_nbasis nn=0;nn<n;nn++) delete [] z[nn];
  delete [] z;
  z=NULL; 
  }
  outf.read((char*)(&N),sizeof(int));
  outf.read((char*)(&n),sizeof(uint_nbasis));
  //preparing space for z
  z=new spsint* [n];
  for (uint_nbasis nn=0;nn<n;nn++) 
  {z[nn]=new spsint [N];
   outf.read((char*)(z[nn]), N*sizeof(spsint));
  }
  //check eof
    if (outf.peek() == EOF) {
        if (outf.eof()) {
            outf.close();
            return 1;
           }
        else {
              outf.close();
        FatalError("Read:ManyBodyStates file-size mismatch ");
                // error
        }        
        }  
 // outf.close();
  return 0;
}



}
//NOTE: This is here because it doesn't need MBS/Map to work so it was "out of place" there. This is probably the best place to put it.
//This converts the many body states from sps_temp to sps. Result is stored back in the same many
//body states object that was given to the function. Assume no sorting is neccessary. (this is due to the
//fact that this will be used when some shells are just omitted so ordering is preserved).
void ConvertMBS(SM::Many_Body_States &st, SM::SingleParticleStates &sps,SM::SingleParticleStates &sps_temp)
{
    std::vector<spsint> key(sps_temp.size());
    spsint tempkey;
    for (int j=0;j<sps_temp.size();j++)
    {
        for (spsint i=0;i<sps.size();i++)
        {
            if (sps_temp[j].N==sps[i].N && sps_temp[j].L==sps[i].L && sps_temp[j].J==sps[i].J && sps_temp[j].Jz==sps[i].Jz && sps_temp[j].Tz==sps[i].Tz) {tempkey=i;break;}
        }
        //            key[i]=SpsIndex(sps_temp,sps,i);
        key[j]=tempkey;
    }
    
    for (uint_nbasis i=0;i<st.n;i++)
    {
        for (int j=0;j<st.N;j++)
            st.z[i][j]=key[st.z[i][j]];
    }
}
//this should be removed later making MBS is different from all these fucntions
#include "MakeMBS.cxx" 
#endif

