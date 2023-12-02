/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file xQN.cxx
\ingroup gp_xsm
*/
#ifndef __XQN__CXX__
#define __XQN__CXX__
#include <sm/QN.cxx>
#include <xml/xmlParser.cxx>
#include <xsm/schema.h>
using xml::XMLNode;
namespace SM{
   int ReadXMLNonAbelianQN(NonAbelianQN &x, XMLNode &xNode){
    x.tag=xNode.getAttribute("name");
	x.J=xNode.getAttribute("j");
	x.T=1;//xNode.getAttribute(xsm::T);
       x.S=1;//this was added by Kostas because MakeNmaxValenceSpace gave different orbitals than reading
	x.N=std::atoi(xNode.getAttribute("n"));
	x.L=std::atoi(xNode.getAttribute("l"));
	//x[i].tag=SPLabel(x[i].N, x[i].L, x[i].J);
	return 0;
}

	std::string WriteXMLNonAbelianQN(XMLNode &xNode, NonAbelianQN &x) {
	addAttribute(xNode,"name",x.tag.c_str());
	addAttribute(xNode,"n",x.N);
    addAttribute(xNode,"j",std::string(x.J).c_str());
    addAttribute(xNode,xsm::T, std::string(x.T).c_str());
    addAttribute(xNode,"l",x.L);
	return xNode.createXMLString(true);
	}
	//this appends QN with whatever is available
int ReadXMLQN(QN &x, XMLNode xNode){ 
	QN xx(x); //save input
    try {getXMLData(x.J,xsm::J,xNode);} catch (...) {x.J=xx.J; }
	try {getXMLData(x.Jz,xsm::Jz,xNode);} catch (...) {x.Jz=xx.Jz;}
	try {getXMLData(x.T,xsm::T,xNode);} catch (...) {x.T=xx.T;}
	try {getXMLData(x.Tz,xsm::Tz,xNode);} catch (...) {x.Tz=xx.Tz;}
	try {getXMLData(x.N,xsm::nucleus_N,xNode);} catch (...) {x.N=xx.N;}
	try {getXMLData(x.L,"L",xNode);} catch (...) {x.L=xx.L;}
	try {std::string tmp; getXMLData(tmp,xsm::P,xNode); x.P=(tmp[0]=='-');}
	catch (...) {x.P=xx.P;}
	try {getXMLData(x.tag,"name",xNode);} catch (...) {x.tag=xx.tag;} 
	/*
	x.J=xNode.getAttribute(xsm::J);
	x.Jz=xNode.getAttribute(xsm::Jz);
	x.T=xNode.getAttribute(xsm::T);
	x.Tz=xNode.getAttribute(xsm::Tz);
	x.N=std::atoi(xNode.getAttribute(xsm::nucleus_N));
	x.L=std::atoi(xNode.getAttribute("L"));
    std::string tmp=xNode.getAttribute(xsm::P);
	x.P=(tmp[0]=='-');
	//x[i].tag=SPLabel(x[i].N, x[i].L, x[i].J);
	*/
	return 0;
}

std::string WriteXMLQN(XMLNode &xNode, QN &x) {
	addAttribute(xNode,xsm::nucleus_N,x.N);
    addAttribute(xNode,xsm::J,std::string(x.J).c_str());
	addAttribute(xNode,xsm::Jz,std::string(x.Jz).c_str());
    addAttribute(xNode,xsm::T, std::string(x.T).c_str());
	addAttribute(xNode,xsm::Tz, std::string(x.Tz).c_str());
    addAttribute(xNode,"L",x.L);
	if (x.P) addAttribute(xNode,xsm::P,"-"); else addAttribute(xNode,xsm::P,"+");
	return xNode.createXMLString(true);
	}

}
#endif //__XQN__CXX__

