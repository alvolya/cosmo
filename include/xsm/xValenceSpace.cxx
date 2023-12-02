/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file xValenceSpace.cxx
\ingroup gp_xsm
*/
#ifndef __XVALENCESPACE__CXX__
#define __XVALENCESPACE__CXX__
#include <xsm/xQN.cxx>
#include <sm/ValenceSpace.cxx>
namespace SM{
std::string WriteXMLValenceSpace(XMLNode &xNode, ValenceSpace &x){
    XMLNode xChildNode;
   for (int i=0;i<x.size();i++) {
   xChildNode=xNode.addChild(xsm::valencespace_orbital);
   WriteXMLNonAbelianQN(xChildNode, x[i]);
   }
	return xNode.createXMLString(true);
}
int ReadXMLValenceSpace(ValenceSpace &x, XMLNode xNode){
   int levels=xNode.nChildNode(xsm::valencespace_orbital);
   if (levels==0) FatalError("There are no s.p. levels");
    x.resize(levels);
	XMLNode xOrbit;
	for (int i=0; i<levels; i++) {
	xOrbit=xNode.getChildNode(xsm::valencespace_orbital,i);
	//this stuff is required to read index correctly
	int ix;
	xml::getXMLData(ix,"index", xOrbit);
	//reading x[ix-1] instead of x[i]
	ReadXMLNonAbelianQN(x[ix-1], xOrbit);
	//x[i].tag=SPLabel(x[i].N, x[i].L, x[i].J);
	}
	return 0;
}
std::string WriteXMLModelSpace(XMLNode &xNode, ModelSpace &x){
  WriteXMLNonAbelianQN(xNode,x);
  return  WriteXMLValenceSpace(xNode,x);
//	return xNode.createXMLString(true);
}

int ReadXMLModelSpace(ModelSpace &x, XMLNode &xNode){
   ReadXMLNonAbelianQN(x,xNode);
   return ReadXMLValenceSpace(x,xNode);
   }

}
#endif //__XVALENCESPACE__CXX__

