/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*!
\author Alexander Volya  <http://www.volya.net>
\file XHH+JJ.cpp
\brief Create Many-body Hamiltonian with JJ-term (xml)
*/
/*SHL M is program to find multipoles
12/24/03 Degenerate model is fine 
02/16/04 particle number is fine
07/10/09 Use new function PPInteraction2sp(VSPX, sps, INT); 
*/


#if __cplusplus <= 199711L
typedef unsigned short spsint;
#else
#include <cstdint>
typedef uint16_t spsint;
#endif
typedef unsigned long long  uint_nbasis;
//#include <math.h>
#define INFO yes
#define STATUS yes
#include <av/statusclass.cxx>
#include <av/argparser.cxx>
#include <xsm/xManyBodyStates.cxx>
#include <xsm/xReadSPInt.cxx>
#include <sm/PPMatrix.cxx>


/*-----------MAIN---------------------------------------*/
int main (int argc, char** argv){
av::argparser args(argc,argv);
args.SetDescription("Create a Hamiltonian matrix and save to file.");
args.SetDefault("XML file containing a parent system description.",2,2);
args.SetFlag("o","output","Output XML file, if provided the default input is unchanged.",0,2);
if(args.ParseArguments()) return 1;
  //unsigned char dfd;
  /*----------Check that interaction file is fine----------*/
  string sysfile=args(1);
 
 
 //make system
XMLNode xMainNode=XMLNode::openFileHelper(sysfile.c_str(), xsm::cosmo); //we start by opening xml input file
 SM::SingleParticleStates sps; Make(sps,xMainNode.getChildNode(xsm::valencespace));
 

 
 XMLNode xSystemNode=xml::SelectNode(xMainNode,xsm::system,"name");
   

 //many-body states
 SM::Many_Body_States st;
SM::Get(st,xMainNode,xSystemNode);
//hamiltonian node
tav::STensor<spsint,double> ESP(2); //output one-body interaction
 tav::STensor<spsint,double> VSPX(4);
 
 //get interaction
xReadSPInt(ESP, VSPX, xSystemNode, xMainNode,sps);


  SM::PPMatrix  HH(st);
  HH.push_back(&ESP);
  HH.push_back(&VSPX);
  //cerr<<VSPX<<endl;
   //XMLNode xNode=xMainNode.addChild("Hamiltonian");
   //SM::WriteXMLPPMatrix(xNode,HH); 
  // cout<<xMainNode.createXMLString(true)<<endl;
 //Pause(); 

//figure out name for output file
string sysname;
 xml::AskAttribute(sysname,"name", "Enter system name: ", xSystemNode);
 string hname=sysname+".HH";
  cerr<<"generating Hamiltonian "<<hname<<endl;
  cerr<<"dimension "<<st.n<<endl;
  HH.Write(hname.c_str()); 

if (args.isSet("o",1)) xMainNode.writeToFile(args("o",1).c_str());
else xMainNode.writeToFile(sysfile.c_str());
 // cout<<xMainNode.createXMLString(true)<<endl;
 return 0;
 }

  







