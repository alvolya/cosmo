/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file Xlevels.cpp
\brief List all levels in a given set of xlm files
*/

//#include <cstdint>
//typedef uint16_t spsint;
//typedef unsigned long long  uint_nbasis;
#include <xsm/xState.cxx>
#include <av/argparser.cxx>

int ListStates(XMLNode xSystemNode) {
int states=xSystemNode.nChildNode(xsm::system_eigenstate);
vector<double> Ediff; //list of theory exp differences of exitation energies
  for (int lst=0;lst<states;lst++) {
	   XMLNode xStateNode=xSystemNode.getChildNode(xsm::system_eigenstate, lst);
	   int id; 
	   double E;
	 //  SM::QN Q;
	  // ReadXMLQN(Q, xStateNode);
	   getXMLData(E,"E",xStateNode);
	   string sE,sEx,sname,sT;
	   getXMLData(sname,"name",xStateNode);
	   getXMLData(sEx,"Ex",xStateNode);
       try {getXMLData(sT,"T",xStateNode); } catch (...) {sT="???";}
	   cout<<sname<<" "<<"T="<<sT<<" "<<sEx<<" "<<E;
     //this is for experimental data
       double Ex;
       getXMLData(Ex,"Ex",xStateNode);
       double Exd; //experimental exitation energy
        //getXMLData(Exd,"Ex",xStateNode.getChildNode("data"));
        try {getXMLData(Exd,"Ex",xStateNode.getChildNode("data"));}
        catch (...) {cout<<endl; continue;}
       cout<<" "<<Exd<<endl;
       Ediff.push_back(Ex-Exd);
	   } //loop over states
    //this prints experimental rms and mean
       if (Ediff.size()!=0) {
       double mean=0.0;
       double sig=0.0;
       for (int s=0;s<Ediff.size();s++) {mean+=Ediff[s]; sig+=Ediff[s]*Ediff[s];}
       mean/=Ediff.size();
       cout<<"Experimental data comparison"<<endl;
       cout<<"Mean (average shift) "<<mean<<endl;
       cout<<"RMS "<<sqrt(sig/Ediff.size()-mean*mean)<<endl;
       }
return 0;
}
/*-----------MAIN---------------------------------------*/
int main (int argc, char** argv){
av::argparser args(argc,argv);
args.SetDescription("List levels and excitation energies in systems/xml files.");
args.SetDefault("XML files containing systems.",2);
args.SetFlag("o","output","Output XML file.",0,2);
if(args.ParseArguments()) return 1;
  //unsigned char dfd;
  /*----------Check that interaction file is fine----------*/
  char sysfile[50];
  char cstmp[256];
  char outftmp[256];
//  if ((argc<2)) {
//    cerr<<"Usage: "<<argv[0]<<" ..xml .... "<<endl;
//	return 1;
//    }
XMLNode xOutMainNodeX=XMLNode::createXMLTopNode("xml",TRUE);
addAttribute(xOutMainNodeX,"version","1.0");
XMLNode xOutMainNodeX1=xOutMainNodeX.addChild("cosmo");
XMLNode xOutMainNode=xOutMainNodeX1.addChild(xsm::nucleus);
//cerr<<args.nArgs("")<<endl;
//return 1;
  for (int argi=1;argi<args.nArgs("");argi++) //loop over all xmlfiles
  {
  XMLNode xMainNode=XMLNode::openFileHelper(args(argi).c_str(), xsm::cosmo); //we start by opening xml input file
   const char *nsys=xsm::nucleus;
  int nsystem=xMainNode.nChildNode(nsys); //number of systems here
   if (nsystem==0) {nsys=xsm::system; nsystem=xMainNode.nChildNode(nsys);} //test systems
  cerr<<"File :"<<argv[argi]<<" systems "<<nsystem<<endl;
  XMLNode xSystemNode;
  for (int nn=0;nn<nsystem;nn++) {
    xSystemNode=xMainNode.getChildNode(nsys, nn);
      int states=xSystemNode.nChildNode(xsm::system_eigenstate);
      if (states==0) continue; //if no states in a given system then continue;
  cerr<<"Model space :"<<xMainNode.getChildNode(nsys, nn).getAttribute("name")<<endl;
  if (Inquire("Include this space ",true)) {
  for (int lst=0;lst<states;lst++) {
	   XMLNode xStateNode=xSystemNode.getChildNode(xsm::system_eigenstate, lst);
	   XMLNode xWFNode=xStateNode.getChildNode(xsm::system_eigenstate_wavefunction);
	   string fname;
	   int id; 
	   double E;
	   SM::QN Q;
	   ReadXMLQN(Q, xStateNode);
	   ReadXMLQN(Q, xWFNode);
	   try {getXMLData(fname,"file",xWFNode); } catch (...) {fname="unknown";}
	   try {getXMLData(id,"id",xWFNode);} catch (...) {id=lst;}
	   getXMLData(E,"E",xStateNode);
	    //SM::AddState(xOutMainNode,Q, E, fname.c_str(), id); //add state to the present system
	   xOutMainNode.addChild(xStateNode.deepCopy()); //this will copy

      
       } //loop over states
  //xOutMainNode.addChild(xMainNode.getChildNode(xsm::system, nn).deepCopy()); //this will copy
  } //ask what system to include
  }//end loop over systems
  }
  SM::OrganizeStates(xOutMainNode);
  ListStates(xOutMainNode);
  //cout<<xOutMainNode.createXMLString(true);
  //xOutMainNode.writeToFile("tmp.xml");

if (args.isSet("o",1)) {
xOutMainNodeX.writeToFile(args("o",1).c_str());
}

 
 return 0;
}









