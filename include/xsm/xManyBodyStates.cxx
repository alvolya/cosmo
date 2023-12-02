/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file xManyBodyStates.cxx
\ingroup gp_xsm
*/
#ifndef __XMANYBODYSTATES__CXX__
#define __XMANYBODYSTATES__CXX__
#include <xsm/xSingleParticleStates.cxx>
#include <sm/ManyBodyStates.cxx>
namespace SM {
/*! \brief Make Many-Body states and save into file, fast method with xml
This is a fast method because particle distribution search is terminated subset has no possibility to satisfy abelian QN
\ingroup gp_mbs
*/
uint_nbasis MakeMBS2(
//int N, int Jz, int Tz, int P, ///<system 
XMLNode xMNode, // system node
SingleParticleStates &sps ///<single-particle states
)
{

  if (xMNode.isEmpty()) FatalError("MakeMBS2: input xml is empty");
  //name many=body states
  char systemname[256]; 
  getXMLData(systemname,"name",xMNode);
  char filename[256]; //name of the output file 
  snprintf (filename, 256, "%s.mbs", systemname); //this is name for mbs
    int N;
	spin Jz, Tz;
  bool P;
   //XMLNode xPNode=xMNode.getChildNode("parameters");
 //if (xPNode.isEmpty()) xPNode=xMNode.addChild("parameters");
    //AskAttribute(x.tag,"name", "Enter name: ", xNode);
	AskAttribute(N,xsm::valence_A, "Enter valence particle number A: ", xMNode);
	Jz=SM::spin(N%2);
	InquireAttribute(Jz,xsm::Jz, "Enter spin projection Jz: ", xMNode);
	Tz=SM::spin(N%2);
	int Nv,Zv;
	try {
	getXMLData(Nv,xsm::valence_N,xMNode);
	getXMLData(Zv,xsm::valence_Z,xMNode);
	Tz=Nv-Zv;
	}
	catch (...) {
	//if we use old xml then this data should be avalable
		InquireAttribute(Tz,xsm::Tz, "Enter valence projection Tz: ", xMNode);
	Zv=(N-Tz)/2;
    Nv=(N+Tz)/2;
	}
	cerr<<"Valence particles N="<<Nv<<" Z="<<Zv<<endl;

	
    char parityz[2];
   AskAttribute(parityz,xsm::P, "Enter parity (+ or -): ", xMNode);
   P=(parityz[0]=='-');
//now we create many-body states node
 XMLNode xNode=xMNode.getChildNode(xsm::system_basisstates);
if (!xNode.isEmpty()) {
  //mbs exist from previous run
  const char *t=xNode.getAttribute("name");
  //if (t!=NULL) {FatalError("file with many-body states is already present")}
  std::cerr<<"basis already present "<<std::endl;
  return 0;
  }
  else xNode=xMNode.addChild(xsm::system_basisstates);
uint_nbasis n;
const char *tox=xMNode.getAttribute(xsm::Nmax); //try to get parameter
if (tox==NULL) {

n = MakeMBS2( filename, N, Jz, Tz, P, sps ); //attribute is not present
}
else {
int Nmax;
getXMLData(Nmax,xsm::Nmax,xMNode);
n = MakeMBS2( filename, N, Jz, Tz, P, Nmax, sps);
}
  
  //XMLNode xFNode=xNode.addChild("file");
  addAttribute(xNode, "name",filename);
  addAttribute(xNode, "size",n);
  addAttribute(xNode, "creator","MakeMBS2");
  
  return n;
}

/*! \brief Make Many-Body states and save into file, fast method.
This is a fast method because particle distribution search is terminated subset has no possibility to satisfy abelian QN
\ingroup gp_mbs
*/
uint_nbasis Make(
SM::Many_Body_States &st,
XMLNode xMainNode, // this main node
XMLNode xSystemNode // this main node  
) {
SM::SingleParticleStates sps; Make(sps,xMainNode.getChildNode(xsm::valencespace));
//XMLNode xSystemNode=xMainNode.getChildNode(xsm::system, -1);
//XMLNode xSystemNode=xml::SelectNode(xMainNode,xsm::system,"name"); //need to select otherwise we always do last one
MakeMBS2(xSystemNode,sps);
//cerr<<xMainNode.createXMLString(true)<<endl;


const char *t=xSystemNode.getChildNode(xsm::system_basisstates).getAttribute("name");
if (t==NULL) FatalError("MakeMBS, cannot find file name in XML")
st.Read(t);
if (st.n==0) {cerr<<"There are no states to work with"<<endl; return 0;}
 cout<<"states "<<st.n<<endl;
 return st.n;
}
    
uint_nbasis Make(
SM::Many_Body_States &st,
XMLNode xMainNode // this main node 
)
{
XMLNode xSystemNode=xml::SelectNode(xMainNode,xsm::system,"name");
return Make(st,xMainNode,xSystemNode);
}

/*! \brief Read or Make Many-Body states and save into file, fast method.
This is a fast method because particle distribution search is terminated subset has no possibility to satisfy abelian QN
\ingroup gp_mbs
*/
uint_nbasis Get(
SM::Many_Body_States &st,
XMLNode xMainNode, // this main node 
XMLNode xSystemNode // this system node 
) {

try{ 
 string cname;
  XMLNode xBasisStates=xml::SelectNode(xSystemNode,xsm::system_basisstates,"name");
  getXMLData(cname,"name", xBasisStates);
 // getXMLData(cname,"name", xSystemNode.getChildNode(xsm::system_basisstates));
 st.Read(cname.c_str());
 }
 catch (...) {SM::Make(st,xMainNode,xSystemNode);}
return st.n;
}
} //end sm namespace
#endif //__XMANYBODYSTATES__CXX__

