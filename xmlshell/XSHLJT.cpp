/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file SHLJT.cpp
\brief Code for analyzing spins and isospins of states, produce *.ev file
*/

#if __cplusplus <= 199711L
typedef unsigned short spsint; //this part is for old c++98
#else
#include <cstdint>
typedef uint16_t spsint;
#endif
typedef unsigned long long  uint_nbasis;
//typedef unsigned long long  uint_nbasis;
//#include <math.h>
//#include <av/sym_sparse.h>
//#include <av/ThreeJSymbol.h>
//#include <sm/ModelSpace.cxx>
#include <av/argparser.cxx>
#define INFO yes
#include <xsm/XSHLJT.cxx>
           





/*-----------MAIN---------------------------------------*/
int main (int argc, char** argv){
av::argparser args(argc,argv);
args.SetDescription("Determine spins and isospins of states and populate xml datables");
args.SetDefault("XML file containing system description.",2,2);
args.SetFlag("o","output","Output XML file, if provided the default input is unchanged.",0,2);
args.SetFlag("n","number","Number of states to put in xml from .EE file ",0,2);
if(args.ParseArguments()) return 1;
//  string sysfile=args(1);
  XMLNode xMainNode=XMLNode::openFileHelper(args(1).c_str(), xsm::cosmo); //we start by opening xml input file
   XMLNode xSystemNode=xml::SelectNode(xMainNode,xsm::system,"name");
 unsigned lstates=10;
  if (args.isSet("n",1)) lstates=args.intarg("n",1);
     xsm::XSHLJT(xSystemNode,xMainNode, lstates);
 //cout<<xMainNode.createXMLString(true)<<endl;

if (args.isSet("o",1)) xMainNode.writeToFile(args("o",1).c_str());
else xMainNode.writeToFile(args(1).c_str());
 return 0;
}









