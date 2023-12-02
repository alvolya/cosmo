/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file xint.cxx
*/
#ifndef __XSPHERICALINT__CXX__
#define __XSPHERICALINT__CXX__
#include <sm/sphericalint.cxx>
namespace SM{
int ReadXMLINT(
      tav::STensor<int,double> &ENT, ///<output: vector data e[1]
      tav::STensor<int,double> &INT, ///<output: sparse-tensor with data V[1][2][3][4][L][T] 
	  XMLNode xMainNode ///< input node that points to xsm::hamiltonian
	  //XMLNode xSystemNode ///< input node that points to system (only get number of particles"
                      ) 
{
	int N;
   AskAttribute(N,xsm::valence_A,"Enter valence particle number: ",xMainNode.getChildNode(xsm::system));
   XMLNode xNode=xMainNode.getChildNode(xsm::hamiltonian) ;
	cerr<<"reading files "<<endl;
	if (xNode.isEmpty()) FatalError("ReadSMinINT:  empty node");
   int files=xNode.nChildNode("file");
   if (files==0) CriticalError("ReadSMinINT(...xml) no files found in xml")
   XMLNode xtmp;
   
	for (int i=0; i<files; i++) {
	xtmp=xNode.getChildNode("file",i);
	//get scaling for this file
	double myscaling[10];
	for (int ii=0;ii<10;ii++) myscaling[ii]=1.0;
	int nscale=xtmp.nChildNode("scale");
	for (int s=0; s<nscale; s++) {
	int rank;
	getXMLData(rank,"rank", xtmp.getChildNode("scale",s));
	double scalerank=atof(xtmp.getChildNode("scale",s).getText());
	cerr<<"Scaling rank= "<<rank<<" with "<<scalerank<<endl;
	myscaling[rank]*= scalerank;
	} 
	
    ReadOxbashINT(ENT,INT,xtmp.getChildNode("name").getText(),N,myscaling[1],myscaling[2]);
	}
	return 0;
}


  int WriteXMLSminINT(
      XMLNode &xNode,
	  std::ifstream &DataStream
                      ) 
{
   int files;
   /*
   int N;
   if (ReadTextData(N,xsm::nucleus_N,DataStream)) {
   cerr<<"Particle number is not defined, enter N= ";
   cin>>N;
   }
*/
   if (ReadTextData(files,"files",DataStream)) CriticalError("There are no interaction files to read");
   char fname[128];
   double spescaling,vscaling;
   XMLNode xChildNode;
   for (int i=0;i<files;i++) {
   DataStream>>fname;
   xChildNode=xNode.addChild("file");
   DataStream>>spescaling;
   DataStream>>vscaling;
   addAttribute(xChildNode,"name",fname);
   addAttribute(xChildNode,"spescaling",spescaling);
addAttribute(xChildNode,"vscaling",vscaling);
   //INF(cerr<<"Reading int file "<<fname<<" "<<spescaling<<" "<<vscaling<<"\n"); 
   //ReadOxbashINT(ENT,INT,fname,N,spescaling,vscaling);
   }
	return 0;
}



}
#endif //__XINT__CXX__

