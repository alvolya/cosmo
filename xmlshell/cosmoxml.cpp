/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file cosmoxml.cpp
\brief Create an xml system description. 
*/
//this revision number will be put into system 
//#ifndef COSMO_REVISION
//#define COSMO_REVISION "xxx"
//#endif
//#define INFO yes;
#include <iostream>
//#include <debug.h>
#include <inquire.cxx>
using tav::Inquire;
#include <av/argparser.cxx>
#include <phy/nuclearnames.cxx>
#include <xsm/cosmoxml.cxx>


#include <av/iofile.cxx> //file with path
//#include <av/wgetfile.cxx> //file with path
using xml::XMLNode;
using namespace std;



int main (int argc, char** argv) {
av::argparser args(argc,argv);
args.SetDescription("XML system creator.");
args.SetDefault("XML file containing xml interaction database or any xml containing system, default: int.xml",0,2);
args.SetFlag("nn","nuclearname","Name of the atomic nucleus e.g. \"18O\" .",0,2);
args.SetFlag("o","output","Output XML file, if not provided user will be asked, if default xml contains only a single interaction the system node will be added.",0,2);
if(args.ParseArguments()) return 1;

string sysfile;
if (args.isSet("",1)) sysfile=args(1);
else sysfile="int.xml";
    XMLNode xMainNode=XMLNode::openFileHelper(av::FileName(sysfile,"SMINTPATH",".xml").c_str(),xsm::cosmo);
    std::cerr<<"Reading int database from \n"<<av::FileName(sysfile,"SMINTPATH",".xml")<<std::endl;
int N,Z;
char nucname[20];
//if (argc>2) std::strcpy(nucname,argv[2]);
//else 
if (args.isSet("nn",1)) std::strcpy(nucname,args("nn",1).c_str());
else AskUser("Enter nucleus name: ", nucname);
 try {
 phy::NuclearNZ(N,Z,nucname);
 }
 catch (...) {
 cerr<<"Name is not understood\n";
 AskUser("Enter neutron number ",N);
 AskUser("Enter proton number ",Z);
 }
 cerr<<"N="<<N<<" Z="<<Z<<endl;
 
 
 //select model space now 
 XMLNode xNodeSP;
 xsm::SelectModelSpace(xNodeSP,xMainNode,N,Z);
//NMIN=  atoi(xNode2.getAttribute(xsm::nucleus_N));
//ZMIN=atoi(xNode2.getAttribute(xsm::nucleus_Z));

XMLNode xNodeINT;
 xsm::SelectModelHamiltonian(xNodeINT, xNodeSP);


//--------------------->start working on output<-------------------------------------

/* output */
int nsp=xMainNode.nChildNode("modelspace"); //number of model spaces
int nint=xNodeSP.nChildNode(xsm::hamiltonian); //number of hamiltonians
XMLNode xOutMainNode;
XMLNode  xOutMainNodeX;
    if ((nsp==0)&&(nint==1)) {
        xOutMainNode=xMainNode;
        xml::DeleteXMLNodes(xOutMainNode,xsm::system);
    }
else
{
//create new
xsm::CoSMoXMLHeader(xOutMainNodeX); //basic header
xOutMainNode=xOutMainNodeX.getChildNode(xsm::cosmo);
xOutMainNode.addChild(xNodeSP.getChildNode(xsm::valencespace)); //add interactions
xOutMainNode.addChild(xNodeINT); //add 
}
//now add system
string nucintname=xsm::AddXMLSystem(xOutMainNode,N,Z);
 
 if (args.isSet("o",1)) { xOutMainNode.writeToFile(args("o",1).c_str()); return 0;}
 //overwrite system 
 if (args.isSet("",1)) { if ((nsp==0)&&(nint==1)) xOutMainNode.writeToFile(args(1).c_str()); return 0;}
  XMLNode xSystem=xOutMainNode.getChildNode(xsm::system);
     string parity;
     parity=xSystem.getAttribute(xsm::P);
  string tmpname=nucintname+parity+".xml";
 //string tmpname=nucintname+".xml";
 Inquire(tmpname, "Name output xml file as ");
 xOutMainNode.writeToFile(tmpname.c_str());
return 0; 
}

