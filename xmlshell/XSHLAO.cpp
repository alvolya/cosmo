/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/** 
\author Alexander Volya  <http://www.volya.net>
\file XSHLAO.cpp
\brief Compute occupation number using Den1B() function.
*/
#if __cplusplus <= 199711L
typedef unsigned short spsint; //this part is for old c++98
#else
#include <cstdint>
typedef uint16_t spsint;
#endif
typedef unsigned long long  uint_nbasis;
#include <av/argparser.cxx>
#include <xsm/xManyBodyStates.cxx>
#include <xsm/xState.cxx>
#include <sm/occupation.cxx>
using namespace SM;
int main(int argc, char** argv) {
av::argparser args(argc,argv);
args.SetDescription("Analyze occupation numbers");
args.SetDefault("XML file containing a parent system description.",2,2);
args.SetFlag("o","output","Output file with plain text",0,2);
if(args.ParseArguments()) return 1;

  string sysfile=args(1);
 XMLNode xMainNode=XMLNode::openFileHelper(sysfile.c_str(), xsm::cosmo); //we start by opening xml input file
 
  XMLNode xValence=xMainNode.getChildNode(xsm::valencespace);
  int levels=xValence.nChildNode(xsm::valencespace_orbital);
  //cerr<<levels<<endl;
 SM::SingleParticleStates sps; Make(sps,xValence);
 
XMLNode xSystemNode=xml::SelectNode(xMainNode,xsm::system,"name");
 if (xSystemNode.isEmpty()) {FatalError("system must exit"); }
 char sysname[256];
 xml::AskAttribute(sysname,"name", "Matrix file name: ", xSystemNode);
 cerr<<" system "<<sysname<<endl;
 //ReadXMLQN(x,xSystemNode);
 //many-body states
 SM::Many_Body_States st;
 SM::Get(st,xMainNode,xSystemNode);
 
    //two systems are read now computing
   // tav::Tensor<2,double> occ(2,x.size());
    vector<double> otx(sps.size());
//std::ofstream ofile;
//ofile.open(argv[2]);
int states=xSystemNode.nChildNode(xsm::system_eigenstate);
if (states==0) FatalError("no states in xml");
cerr<<" states "<<states<<endl;
XMLNode xStateNode;
double E;
tav::Tensor<1,double> z;
       for (int lst=0;lst<states;lst++) {
	   xStateNode=xSystemNode.getChildNode(xsm::system_eigenstate, lst);
	   cout<<xStateNode.createXMLString(true)<<endl;
	   int id; 
	   char fname[256];
	   XMLNode xWFNode=xStateNode.getChildNode(xsm::system_eigenstate_wavefunction);
	   getXMLData(fname,"file",xWFNode);
	   getXMLData(id,"id",xWFNode);
	   //getXMLData(E,xsm::E,xStateNode);
	   EnergyStateRead(E, z, fname, id);
	  // cerr<<fname<<endl;
    //   cout<<E<<endl;
       //this is overlap
    

SM::occupation(otx,z.elem,st);
 XMLNode xOrbital;
 for (int orb=0;orb<levels;orb++) {
 xOrbital=xValence.getChildNode(xsm::valencespace_orbital, orb);
 double occN=0.0;
 double occZ=0.0;
 int inx;
 xml::getXMLData(inx, "index", xOrbital);
 inx--;
 for (int i=0;i<sps.size();i++) if (sps[i].level== inx) if (sps[i].Tz==1) occN+=otx[i]; else occZ+=otx[i];
 XMLNode xDen=xStateNode.addChild("occupation");
 //if (!xDen.isEmpty()) FatalError("occupation numbers are already present");
// xDen=xStateNode.addChild("occupation");
 addAttribute(xDen,"name", xOrbital.getAttribute("name"));
 addAttribute(xDen,xsm::nucleus_N,occN);
 addAttribute(xDen,xsm::nucleus_Z,occZ);
 }

       }
	 //  xMainNode.writeToFile(sysfile);
if (args.isSet("o",1)) xMainNode.writeToFile(args("o",1).c_str());
else xMainNode.writeToFile(args(1).c_str());
	//   cout<<xMainNode.createXMLString(true)<<endl;
  return 0;
}

