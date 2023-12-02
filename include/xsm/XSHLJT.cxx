/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file XSHLJT.cxx
\brief Create xml database of states
\ingroup gp_xsm
*/
#ifndef __XSHLJT__CXX__
#define __XSHLJT__CXX__


#include <algorithm>
#include <sm/SingleParticleStates.cxx>
#include <xsm/xManyBodyStates.cxx>
//will include system and single particle states
 
//#include <sm/JJVTimesVector.cxx>
#include <xsm/xState.cxx>
#include <sm/parityclass.cxx>
using SM::spin;
//#include <sm/supplemental.cxx>
//#include "MakeOX.cxx"
//#include <sm/MakeSPO.cxx>
//#include "QMatrix.cxx"
//#include "QQMatrix.cxx"
#include <sm/spherical2sp.cxx>
#include <sm/PPMatrix.cxx>  
           
using SM::State;
#define INFO yes
namespace xsm {
int XSHLJT(
XMLNode xSystemNode, ///< [i/o] system node
XMLNode xMainNode, ///< [i] main node for sps
unsigned lstates ///< number of states to process 
)
{
SM::SingleParticleStates sps; Make(sps,xMainNode.getChildNode(xsm::valencespace));
 
 //XMLNode xSystemNode=xMainNode.getChildNode(xsm::system,-1); //get last child
 
 XMLNode xMatrixNode=xSystemNode.getChildNode(xsm::system_matrix,-1); //get last child
 if (xMatrixNode.isEmpty()) FatalError("You need Hamiltonian matrix");
 string sysname;
 //-------if system defines T or J then we use them
 spin xJ;
 try {xml::getXMLData(xJ, "J", xSystemNode); } catch (...) {xJ=-1;}
 if (xJ!=-1) cerr<<"States only with J="<<xJ<<" are processed \n";
 spin xT;
 try {xml::getXMLData(xT, "T", xSystemNode); } catch (...) {xT=-1;}
 if (xT!=-1) cerr<<"States only with T="<<xT<<" are processed \n";
 
 xml::AskAttribute(sysname,"name", "Matrix file name: ", xMatrixNode);
  double jscale=0.0;
 try { xml::getXMLData(jscale, "JJ", xMatrixNode); } catch (...) {jscale=0.0;} 
 cerr<<"J+J- scaling is "<<jscale<<endl;
   double tscale=0.0;
 try { xml::getXMLData(tscale, "TT", xMatrixNode); } catch (...) {tscale=0.0;}
 cerr<<"T+T- scaling is "<<tscale<<endl;
 string dataformat;
 if (xml::ReadXMLData(dataformat,"format",xMatrixNode)) {CriticalError("matrix is not in pp format");}
 else if (dataformat!="pp") {CriticalError("matrix is not in pp format");}
 
 SM::QN x; //unfinished
 ReadXMLQN(x,xSystemNode);
 //many-body states
 SM::Many_Body_States st;
SM::Get(st,xMainNode,xSystemNode);
 
 tav::STensor<spsint,double> JO(2);
   SM::SPJPlus(JO, sps);
 tav::STensor<spsint,double> TO(2);
   SM::SPTPlus(TO, sps);
   SM::PPMatrix JJ(st); JJ.push_back(&JO);
   SM::PPMatrix TT(st); TT.push_back(&TO);
   

string cstmp=sysname+".EE";
ifstream outf;


 outf.open(cstmp.c_str(), ios::binary);
 cout<<"reading file "<<cstmp<<endl;
 unsigned lowlyingstates=1;
 rw::Read(lowlyingstates, outf);
 cerr<<"states found "<<lowlyingstates<<endl;

 tav::Tensor<1,double> z;
 if ((lstates>0)&&(lstates<lowlyingstates)) lowlyingstates=lstates;
 else
 if (lowlyingstates>20) {cerr<<"There are "<<lowlyingstates<<", how many do you want to process? "; 
     int itmp;
     cin>>itmp;
     if (itmp>0) lowlyingstates=itmp; 
 }
 cerr<<"processing states "<<lowlyingstates<<endl;
 for (int lst=0;lst<lowlyingstates;lst++) {
 SM::QN Q(x); //quantum numbers of my state
    //------------------------------------------------------
   double E;
   double deltaE;
   double deltaV;
   rw::Read(E, outf); //read energy
   rw::Read(deltaE, outf); //read energy
   rw::Read(deltaV, outf); //read energy
   z.Read(outf);
  double je=0.;
  je=0.0;
  je=JJ.Sqr(z); 
  cerr<<"State: "<<lst<<" E= "<<E;
  cerr<<" JJ= "<<je<<endl;
    SM::spin J;
    J= SM::JJSpin(je,x.Jz);
    if ((xJ!=-1)&&(xJ!=J)) continue;
 if (fabs(jscale)>1E-8) {
    E-=jscale*((J-x.Jz)*(J+x.Jz+2)/4.0);
    cerr<<"Applying shift js="<<jscale<<" s="<<jscale*((J-x.Jz)*(J+x.Jz+2)/4.0)<<endl;
    }
   Q.J=J;
   je=TT.Sqr(z); 
   cerr<<" TT= "<<je<<" "<<endl;
   J= SM::JJSpin(je,x.Tz);
       if ((xT!=-1)&&(xT!=J)) continue;
    if (fabs(tscale)>1E-8) {
    E-=tscale*((J-x.Tz)*(J+x.Tz+2)/4.0);
    cerr<<"Applying shift ts="<<tscale<<" s="<<tscale*((J-x.Tz)*(J+x.Tz+2)/4.0)<<endl;
    }
   Q.T=J;
   //SM::AddState(xMainNode.getChildNode(xsm::system,0),Q, E, cstmp, lst); 
   XMLNode xStateNode= SM::AddState(xSystemNode,Q, E, cstmp.c_str(), lst); //add state to the present system
    updateAttribute(xStateNode,"dE",deltaE);
 }
// SM::OrganizeStates(xMainNode.getChildNode(xsm::system,0));
 SM::OrganizeStates(xSystemNode);
 return 0;
}

}//end namespace xsm









#endif //__XSHLJT__CXX__

