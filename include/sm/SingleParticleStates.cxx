/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file SingleParticleStates.cxx
\brief M-scheme s.p. states 
*/
#ifndef __SPS__CXX__
#define __SPS__CXX__
#include <sm/ValenceSpace.cxx>
#include <sm/parityclass.cxx>
#include <vector>
#include <map>
#include <limits>  //check spsint
#include <algorithm>    // std::sort
using std::vector;
using std::map; //use to quickly count levels
/*
**\typedef spsint 
\brief spsint is type, sufficient to record all single-particle locations 
  some codes do not define spsint, it is approprate to place it here. 
  \todo if this is defined in SingleParticleStates.cxx then we can not overwrite this, need compiler instruction?
present definition must be global because of cioperators to work
different isospins are not possible
 */
//typedef unsigned char spsint; 
 
namespace SM{    
/*! Data structure for spherical s.p. state
*/
  struct SPS: QN{
      //!particle number
    int level; ///< reference to a s.p. level
  };

std::ostream& operator << (std::ostream& stream, const SPS q)
{
  stream<<q.N<<'\t';
  stream<<q.J<<'\t';
  stream<<q.Jz<<'\t';
  stream<<q.T<<'\t';
  stream<<q.Tz<<'\t';
  stream<<q.P<<'\t';
  stream<<q.level<<'\n';
        return stream;
}

template <class cnumber>  
std::istream& operator >> (std::istream& stream, SPS &q)
{
  stream>>q.N;
  stream>>q.J;
  stream>>q.Jz;
  stream>>q.T;
  stream>>q.Tz;
  stream>>q.P;
  stream>>q.level;
        return stream;
}

/*! \brief Returns the number of quanta in a single particle state*/
    int Quanta(SPS x)
    {
        return x.N*2 + x.L;
    }

/*! Array of m-scheme s.p. states
\typedef SingleParticleStates
\brief Vector of s.p. states 
 
This is an array of SPS, each element corresponds to an m-scheme 
position, and refers to all quantum numbers related to this position
*/
typedef std::vector<SM::SPS> SingleParticleStates;



 /*! 
\brief Function that creates Array of s.p. states from valence space
This order is special doing neutrons first and then protons. Acceptible for rejection technology. 
*/
 int Make(
 SingleParticleStates &sps, ///<[output] S.P. states
 ValenceSpace xq ///< quantum numbers 
 ) {
int levels=xq.size();
int n=0;
int itmp;
for (itmp=0; itmp<levels; itmp++) n+=(xq[itmp].J+1)*(xq[itmp].T+1);
  if (std::numeric_limits<spsint>::max()<n) {
   FatalError("spsint variable type is too small !!!!\n The "<<n<<" single particle states is more than \n"<<std::numeric_limits<spsint>::max()<<" single particle states: "<<(unsigned int)(n)<<"\n");}
    //else {cerr<<"Maximum value for <spsint> "<<std::numeric_limits<spsint>::max()<<endl;}
 
   sps.resize(n); //make of an appropriate size
itmp=0;
//for (int tt=1;tt>=-1;tt=tt-2) 
//(int tt=-x.t[lt];tt<=x.t[lt];tt=tt+2)
for (int tt=xq[0].T ;tt>=-xq[0].T;tt=tt-2)
for (int lt=0; lt<levels; lt++) 
   for (int m=-xq[lt].J;m<=xq[lt].J;m=m+2) 
{
sps[itmp].J=xq[lt].J;
sps[itmp].T=xq[lt].T;
sps[itmp].L=xq[lt].L;
sps[itmp].S=xq[lt].S;
sps[itmp].N=xq[lt].N;
sps[itmp].P=(xq[lt].L&1);  //orbital l
sps[itmp].Tz=tt;
sps[itmp].Jz=m;
 sps[itmp].level=lt;
itmp++;
 }
 // cout<<"single particle states: "<<(unsigned int)(n)<<endl;
 return n;
}



/*!Function to find a give s.p. state using level, Jz and Tz
\return number corresponding to s.p. m-scheme state 
*/
spsint Find(
SingleParticleStates &sps,///<reference to an array of states 
int level, ///< s.p. level of interest
int Jz, ///< s.p. Jz projection
int Tz ///<s.p. Tz projection
) {
  spsint itmp=0;
  for (itmp=0;itmp<sps.size();itmp++) {
    if (sps[itmp].level!=level) //continue;
     {itmp+=spsint(sps[itmp].J); continue; /*add 2J moving toward end and continue */}
    if (sps[itmp].Jz!=Jz) continue;
    if (sps[itmp].Tz!=Tz) continue;
    return itmp;
}
return itmp;
}

int CountLevels(SM::SingleParticleStates &sps) {
	std::map<int, int> tmpmap;
	for (int is=0;is<sps.size();is++) tmpmap[sps[is].level]=0; //identify all levels
	return tmpmap.size();
	}

/*
 \fn CoreQuanta
 \brief Count number of core quanta (minimum) for a given numbar of particles A
 quanta can be separated for protons and neutrons
 */
int CoreQuanta(int A, ///<Number of particles
               int T, ///< Isospin, >0 neutrons, <0, protons =0 both
               SM::SingleParticleStates &sps) {
std::vector<int> quanta;
    
    if (T==0) for (int is=0;is<sps.size();is++) quanta.push_back(Quanta(sps[is]));
    if (T>0) for (int is=0;is<sps.size();is++) if (sps[is].Tz>0) quanta.push_back(Quanta(sps[is]));
    if (T<0) for (int is=0;is<sps.size();is++) if (sps[is].Tz<0) quanta.push_back(Quanta(sps[is]));
    if (A>quanta.size()) return -1; //return error
    std::sort (quanta.begin(), quanta.begin()); //sort in ascending order
    int QQ=0;
    for (int i=0;i<A;++i) QQ+=quanta[i];
return QQ;
}


/*! \fn StateAQN
 \brief Obtain Abelian QN from a given pp-state (array of s.p. indexes)
 */
    AbelianQN StateAQN(SM::SingleParticleStates &sps, unsigned n, spsint *mbs) {
        AbelianQN x; x.Jz=0; x.Tz=0; x.P=0;
        for (int i=0;i<n;i++) {x=x+sps[mbs[i]]; }
       // for (int i=0;i<n;i++) {x.Jz=+sps[mbs[i]].Jz; x.Tz+=sps[mbs[i]].Tz; x.P^=sps[mbs[i]].P;}
        return x;
    }
 

/*!\brief Mapping vector from one sps to another set of sps 
If for initial state mapping does not exits 
std::numeric_limits<spsint>::max() is returned, we use spsint(-1) for that.
*/
 class SPSKeyMap: public std::vector<spsint> {
 //constructor is the only funciton here
 public:
   SPSKeyMap(SM::SingleParticleStates& sfinal,SM::SingleParticleStates& sinitial)
        : std::vector<spsint> (sinitial.size(),spsint(-1)) //initialize vector with final size value
   {
    for (auto i = 0; i < sinitial.size(); i++)
        for (auto j = 0; j < sfinal.size(); j++)
            if (sinitial[i] == sfinal[j])
            {
                (*this)[i] = j;
                break;//found no need to continue.
        }
   }
 };
    
} //end of namespace SM

#endif

