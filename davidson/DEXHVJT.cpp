/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*!
\author Alexander Volya  <http://www.volya.net>
\file DEXHVJT.cpp
\brief Davidson diagonalization of Virtual Hamiltonian (generated on fly) 
*/
/*SHL M is program to find multipoles
12/24/03 Degenerate model is fine 
02/16/04 particle number is fine
07/10/09 Use new function PPInteraction2sp(VSPX, sps, INT); 
*/

#if __cplusplus <= 199711L
typedef unsigned short spsint; //this part is for old c++98
#else
#include <cstdint>
typedef uint16_t spsint;
#endif
typedef unsigned long long  uint_nbasis;
//#include <math.h>
#define INFO yes
#include <av/argparser.cxx>
#include <av/statusclass.cxx>
#define CSTATUS
#define FSTATUS
#include <xsm/DEXHV.cxx>
#include <xsm/XSHLJT.cxx>

/*-----------MAIN---------------------------------------*/
int main (int argc, char** argv){
av::argparser args(argc,argv);
args.SetDescription("Davidson or eigen diagonalization of Virtual Hamiltonian (generated on fly).");
args.SetDefault("XML file containing a system description",2,2);
args.SetFlag("i","iterations","Set the maximum number of iterations.",0,2);
args.SetFlag("n","number","Set the desired number of states.",0,2);
args.SetFlag("b","block","Set the block size.",0,2);
args.SetFlag("f","file","File name where initial set of fectors is stored",0,2);
args.SetFlag("o","output","Output XML file, if provided the default input is unchanged.",0,2);
args.SetFlag("l","limit","Limit under which exact diagonalization is used. (default is 1000, set to zero for exact diagonalization regardless of size)",0,2);

if(args.ParseArguments()) return 1;

int iters=0; int bloc=10; int nstate=0;
    int limit = 1000;
    string invectorf="NULL"; //file for in vectors
if (args.isSet("i",1)) iters=args.intarg("i",1);
if (args.isSet("n",1)) nstate=args.intarg("n",1);
if (args.isSet("b",1)) bloc=args.intarg("b",1);
if (args.isSet("f",1)) { invectorf=args("f",1); }
    if (args.isSet("l",1)) limit = args.intarg("l",1);
//cerr<<iters<<" "<<nstate<<" "<<bloc<<endl;
if (bloc<nstate) bloc=nstate; //block must be bigger then number of saved states
  //unsigned char dfd;
  /*----------Check that interaction file is fine----------*/
  string sysfile=args(1);
 
 //make system
    XMLNode xMainNode=XMLNode::openFileHelper(sysfile.c_str(), xsm::cosmo); //we start by opening xml input file
 XMLNode xSystemNode=xml::SelectNode(xMainNode,xsm::system,"name");
   
   //diagonalize the hamiltonian
  xsm::DXHV(xSystemNode, xMainNode, iters, bloc, nstate, invectorf, limit);
  //populate database
  
  xsm::XSHLJT(xSystemNode,xMainNode, nstate);
  
  if (args.isSet("o",1)) xMainNode.writeToFile(args("o",1).c_str());
else xMainNode.writeToFile(sysfile.c_str());
  //xMainNode.writeToFile(sysfile);
 // cout<<xMainNode.createXMLString(true)<<endl;
 return 0;
}

  







