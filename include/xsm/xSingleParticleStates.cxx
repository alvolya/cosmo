/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file xSingleParticleStates.cxx
\ingroup gp_xsm
*/
#ifndef __XSINGLEPARTICLESTATES__CXX__
#define __XSINGLEPARTICLESTATES__CXX__
#include <xsm/xValenceSpace.cxx>
#include <sm/SingleParticleStates.cxx>
#include <map>
#include <limits>  //check spsint
using std::map;
namespace SM{
/*! 
\brief Function that creates Array of s.p. states from valence space
This order is special doing neutrons first and then protons. Acceptible for rejection technology. 
The data in xml is sorted by index
*/
 int Make(
 SingleParticleStates &sps, ///<[output] S.P. states
 XMLNode xNode ///<[input] xml node pointing to valencespace;
 ) {
if (xNode.isEmpty()) CriticalError("Make(sps,xml): empty xml on input");
int levels=xNode.nChildNode(xsm::valencespace_orbital);
if (levels==0) CriticalError("Make(sps,xml): no orbitals defined in xml file");
//cerr<<"got levels "<<levels<<endl;
std::vector<SM::SPS> protons;
std::vector<SM::SPS> neutrons;
//iterate over all orbitals in xml
XMLNode xtmpNode;
//first we need to sort the data in order of inex
std::map<int,int> orderinx; //make array for sorting
for (int ll=0;ll<levels;ll++) {
xtmpNode=xNode.getChildNode(xsm::valencespace_orbital,ll);
int indx;
getXMLData(indx,"index",xtmpNode); //read from xml index for a given level
orderinx[indx]=ll;
}
//map comes to be already sorted so we are good
//we iterate over all levels using map which automatically sorts them
for (std::map<int,int>::iterator linx=orderinx.begin();linx!=orderinx.end();linx++) {
xtmpNode=xNode.getChildNode(xsm::valencespace_orbital,linx->second); //we read the child according to map
SPS x;
 getXMLData(x.tag,"name",xtmpNode);
 //x.tag=xtmpNode.getAttribute("name");
 //x.J=xtmpNode.getAttribute("j");
 getXMLData(x.J,"j",xtmpNode);
 x.N=0; InquireXMLData(x.N,"n",xtmpNode); //read main quantum number, default 0
 x.L=0; InquireXMLData(x.L,"l",xtmpNode); //read orbital quantum number, default 0
 // x.N=std::atoi(xtmpNode.getAttribute("n"));
 //x.L=std::atoi(xtmpNode.getAttribute("l"));
 getXMLData(x.level,"index",xtmpNode);
 x.level--; //fortran conversion
 	
 //x.level=std::atoi(xtmpNode.getAttribute("index"))-1; //note this is fortran to c conversion
 //other standard things
 x.S=1;
 SM::parity PP=(x.L&1);
 InquireXMLData(PP,xsm::P,xtmpNode);
 x.P=bool(PP); //do not use x.P directy because it is of different type (int)
 x.T=1; 
 InquireXMLData(x.T, xsm::T, xtmpNode); //uses default x.T=1
 //Isospin is done, try projections
 std::string type="pn";
 try {getXMLData(x.Tz, xsm::Tz, xtmpNode);
 if (x.Tz==1) type="n";
 if (x.Tz==-1) type="p";
 }
 catch (...) {
 type="pn";
 }//we were not able to determine type 
//assume pn by default set by Tz above
InquireXMLData(type,"type",xtmpNode);
 //type=xtmpNode.getAttribute("type");
 //work on m projection
 for (int m=-x.J;m<=x.J;m=m+2) {
 x.Jz=m;
 if (x.T==0) {x.Tz=0;  neutrons.push_back(x); continue;}
 if (type=="pn") {x.Tz=1; neutrons.push_back(x); x.Tz=-1; protons.push_back(x);}
 if (type=="p") {x.Tz=-1; protons.push_back(x);}
 if (type=="n") {x.Tz=1; neutrons.push_back(x);}
 }
 } //go over levels

//we now have proton and neutron levels
 int n=protons.size()+neutrons.size();
 if (std::numeric_limits<spsint>::max()<n) {
   FatalError("spsint variable type is too small !!!!\n The "<<n<<" single particle states is more than \n"<<std::numeric_limits<spsint>::max()<<" single particle states: "<<(unsigned int)(n)<<"\n");}
    else {cerr<<"Maximum value for <spsint> "<<std::numeric_limits<spsint>::max()<<endl;}
   sps.resize(n); //make of an appropriate size
//incert everything 
  int itmp=0;
  for (int i=0;i<neutrons.size();i++) {sps[itmp]=neutrons[i]; itmp++;}
  for (int i=0;i<protons.size();i++) {sps[itmp]=protons[i]; itmp++;}
 // cout<<"single particle states: "<<(unsigned int)(n)<<endl;
 return n;
}
}

#endif //__XSINGLEPARTICLESTATES__CXX__

