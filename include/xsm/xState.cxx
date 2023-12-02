/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file xState.cxx
\ingroup gp_xsm
*/
#ifndef __XSTATE__CXX__
#define __XSTATE__CXX__
#include <xsm/xQN.cxx>
#include <sm/State.cxx>
#include <map>
using std::map;
namespace SM{
/*Read a wave function for a given state does not need a State class
*/
template <class wfclass>
int XMLReadWaveFunction(
double &Energy, //< [0] energy of the state
wfclass &z, //<[0] the wave function 
XMLNode xWFNode //pointer to wf node
) {
	   std::string fname;
	   int id;
	   if (xWFNode.isEmpty()) CriticalError("No wave funciton found for a state")
	   getXMLData(fname,"file",xWFNode);
	   getXMLData(id,"id",xWFNode);
	   //getXMLData(E,xsm::E,xStateNode);
	   SM::EnergyStateRead(Energy, z, fname.c_str(), id);
	return id;
}

int ReadXMLState(SM::State &PT, xml::XMLNode xStateNode) {
	   SM::ReadXMLQN(PT, xStateNode);
	   XMLNode xWFNode=xStateNode.getChildNode(xsm::system_eigenstate_wavefunction);
	   SM::ReadXMLQN(PT, xWFNode); //additional QN are in WF
	   //getXMLData(E,xsm::E,xStateNode);
	   XMLReadWaveFunction(PT.E,PT.z,xWFNode);
	   getXMLData(PT.tag,"name", xStateNode);
	   return 0;
}



/*Read a wave function for a given state does not need a State class
*/
template <class wfclass>
int XMLEnergyStateRead(
double &Energy, //< [0] energy of the state
wfclass &z, //<[0] the wave function 
XMLNode xSystemNode, //pointer to system
int i//which state to read 
) {
XMLNode xStateNode=xSystemNode.getChildNode(xsm::system_eigenstate, i);
	   //cout<<xStateNode.createXMLString(true)<<endl;
	   int states=xSystemNode.nChildNode(xsm::system_eigenstate);
	   if (i>=states) CriticalError("There is no such state");
	   XMLNode xWFNode=xStateNode.getChildNode(xsm::system_eigenstate_wavefunction);
	   if (xWFNode.isEmpty()) CriticalError("No wave funciton found for a state")
	   XMLReadWaveFunction(Energy,z,xWFNode);
	return states;
}

/*! \brief Search for a many body state search with given parameters
*/
XMLNode SearchState(XMLNode xSystemNode, spin JJ, spin TT, double E, double DE=0.001) {
int nsp=xSystemNode.nChildNode(xsm::system_eigenstate); //number of states
//cerr<<JJ<<endl;
//cerr<<nsp<<endl;
for (int nn=0;nn<nsp;nn++) {
XMLNode xNodeSelect=xSystemNode.getChildNode(xsm::system_eigenstate,nn);
spin myJ;
getXMLData(myJ,xsm::J,xNodeSelect);
//cerr<<myJ<<" J="<<JJ<<endl;
if (myJ!=JJ) continue;
getXMLData(myJ,xsm::T,xNodeSelect);
if (myJ!=TT) continue;
double myE;
getXMLData(myE,xsm::E,xNodeSelect);
if (std::abs(myE-E)>DE) continue;
return xNodeSelect;
}
//cerr<<"Node note found;\n";
CriticalError("State not found");
}


/*! \brief Search state by energy only, states must be ordered
\return integer, positive if found match, negative if not found, counting starts from 1  
*/
int SearchState(XMLNode xSystemNode, double E, double DE=0.001) {
int nsp=xSystemNode.nChildNode(xsm::system_eigenstate); //number of states
//cerr<<JJ<<endl;
//cerr<<nsp<<endl;
int nn=1;
for (nn=1;nn<=nsp;nn++) {
XMLNode xNodeSelect=xSystemNode.getChildNode(xsm::system_eigenstate,nn+1);
double myE;
getXMLData(myE,xsm::E,xNodeSelect);
if (std::abs(myE-E)<DE) return nn; //energy match
if (E<myE) return -nn;
}
//cerr<<"Node note found;\n";
//CriticalError("State not found");
return -nn;
}

XMLNode AddState(XMLNode xSystemNode,SM::QN Q, double E, const char *filename, int id) {
// define XMLNode and check if the state is already in 
XMLNode xStateNode;
try {
xStateNode=SearchState(xSystemNode,Q.J, Q.T, E);
}
//-----------------------------



catch (...) {
//  cerr<<"adding new node\n";
   xStateNode=xSystemNode.addChild(xsm::system_eigenstate);
   updateAttribute(xStateNode,"name",string(Q.J));
   updateAttribute(xStateNode,xsm::J,string(Q.J));
   if (Q.P) updateAttribute(xStateNode,xsm::P,"-"); else updateAttribute(xStateNode,xsm::P,"+");
   updateAttribute(xStateNode,xsm::T,string(Q.T));
    addAttribute(xStateNode,xsm::E,E);
}

//add new wf
//see if such wave function is already present
XMLNode xWFNode=xStateNode.getChildNodeWithAttribute(xsm::system_eigenstate_wavefunction,"file",filename); 
if (!xWFNode.isEmpty()) {
int id2;
getXMLData(id2, "id" ,xWFNode); 
if (id==id2) return xStateNode; //the wf is already in 
}
//actually add new wf

xWFNode=xStateNode.addChild(xsm::system_eigenstate_wavefunction);
  updateAttribute(xWFNode,xsm::Jz,Q.Jz);
  updateAttribute(xWFNode,xsm::Tz,Q.Tz);
  updateAttribute(xWFNode,"file",filename);
  updateAttribute(xWFNode, "id", id);
return xStateNode;
}

/*
Organize and name many-body states in the system
*/
int NameStates(XMLNode xStateNode) {
int nsp=xStateNode.nChildNode(xsm::system_eigenstate); //number of states
std::vector<std::pair<double, int> > energy(nsp);
for (int nn=0;nn<nsp;nn++) {
XMLNode xNodeSelect=xStateNode.getChildNode(xsm::system_eigenstate,nn);
/*
spin myJ;
getXMLData(myJ,xsm::J,xNodeSelect);
//cerr<<myJ<<" J="<<JJ<<endl;
if (myJ!=JJ) continue;
getXMLData(myJ,xsm::T,xNodeSelect);
if (myJ!=TT) continue;
*/
double myE;
getXMLData(myE,xsm::E,xNodeSelect);
energy[nn].first=myE;
energy[nn].second=nn;
}
std::sort(energy.begin(), energy.end());
std::map< std::pair<int,bool> ,int> spinorder;
//name states
for (int in=0;in<nsp;in++) {
int nn=energy[in].second;
XMLNode xNodeSelect=xStateNode.getChildNode(xsm::system_eigenstate,nn); 
spin myJ;
getXMLData(myJ,xsm::J,xNodeSelect);
SM::parity PP;
getXMLData(PP,xsm::P,xNodeSelect);
//bool PP=(parity[0]=='-');
std::pair<int,bool> JP(myJ,PP);
if (spinorder.find(JP)==spinorder.end()) spinorder[JP]=1; else spinorder[JP]++;
char statename[32];
sprintf(statename,"%s%c(%i)",string(myJ).c_str(),char(PP),spinorder[JP]);
//update label and energy of a state
//cerr<<statename<<" "<<energy[in].first-energy[0].first<<endl;
updateAttribute(xNodeSelect,"name",statename);
updateAttribute(xNodeSelect,xsm::Ex,energy[in].first-energy[0].first);
}
return 0;
}

/*
Organize and name many-body states in the system
*/
int OrganizeStates(XMLNode xStateNode) {
if (xStateNode.isEmpty()) {std::cerr<<"OrganizeStates:: empty node \n"; return 0;}
int nsp=xStateNode.nChildNode(xsm::system_eigenstate); //number of states
std::vector<std::pair<double, int> > energy(nsp);
for (int nn=0;nn<nsp;nn++) {
XMLNode xNodeSelect=xStateNode.getChildNode(xsm::system_eigenstate,nn);
/*
spin myJ;
getXMLData(myJ,xsm::J,xNodeSelect);
//cerr<<myJ<<" J="<<JJ<<endl;
if (myJ!=JJ) continue;
getXMLData(myJ,xsm::T,xNodeSelect);
if (myJ!=TT) continue;
*/
double myE;
getXMLData(myE,xsm::E,xNodeSelect);
energy[nn].first=myE;
energy[nn].second=nn;
}
std::sort(energy.begin(), energy.end());
std::map< std::pair<string,string> ,int> spinorder;
//name states
XMLNode xTemporary=xStateNode.addChild("temporary");
for (int in=0;in<nsp;in++) {
int nn=energy[in].second;
XMLNode xNodeSelect=xStateNode.getChildNode(xsm::system_eigenstate,nn); 
string myJ;
try {getXMLData(myJ,xsm::J,xNodeSelect);} catch (...) {myJ="??";}
string PP;
try {getXMLData(PP,xsm::P,xNodeSelect); } catch (...) {PP="?";}
//bool PP=(parity[0]=='-');
std::pair<string,string> JP(myJ,PP);
if (spinorder.find(JP)==spinorder.end()) spinorder[JP]=1; else spinorder[JP]++;
char statename[32];
sprintf(statename,"%s%c(%i)",string(myJ).c_str(),char(PP[0]),spinorder[JP]);
//update label and energy of a state
//cerr<<statename<<" "<<energy[in].first-energy[0].first<<endl;
updateAttribute(xNodeSelect,"name",statename);
updateAttribute(xNodeSelect,xsm::Ex,energy[in].first-energy[0].first);
xTemporary.addChild(xNodeSelect.deepCopy());
}
//delete all states
//cerr<<xStateNode.createXMLString(true)<<endl;
//cerr<<xTemporary.createXMLString(true)<<endl;
for (int in=0;in<nsp;in++) (xStateNode.getChildNode(xsm::system_eigenstate)).deleteNodeContent(); 
for (int in=0;in<nsp;in++) xStateNode.addChild(xTemporary.getChildNode(xsm::system_eigenstate,in).deepCopy()); //move from xTemporary
xTemporary.deleteNodeContent();
//cerr<<xStateNode.createXMLString(true)<<endl;
return nsp;
}

}
#endif //__XSTATE__CXX__

