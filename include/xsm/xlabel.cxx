/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file xlabel.cxx
\ingroup gp_xsm
*/
#ifndef __XLABEL__CXX__
#define __XLABEL__CXX__

#include <xml/xmlParser.cxx>
using xml::XMLNode;
#include <sm/spinclass.cxx> //contains string stream
using std::stringstream;
using SM::spin;
namespace SM{
std::string LabelValenceSpace(xml::XMLNode xModelSpace) {
  int levels=xModelSpace.nChildNode(xsm::valencespace_orbital);
   if (levels==0) FatalError("There are no s.p. levels");
	spin J;
	std::stringstream retstr;
	for (int i=0; i<levels; i++) {
	xml::getXMLData(J,"j", xModelSpace.getChildNode(xsm::valencespace_orbital,i));
	if (i==0) {if (int(J)&1) retstr<<"f"; else retstr<<"b";}
	retstr<<"l"<<int(J)+1;
	}
	return retstr.str();
}

std::string LabelSystem(xml::XMLNode xSystem) {
string readstring;
std::stringstream returnstring;
if (!ReadXMLData(readstring,xsm::nucleus_N,xSystem)) {returnstring<<"N"<<readstring;}
spin J;
if (!ReadXMLData(J,xsm::Jz,xSystem)) {returnstring<<"M"<<int(J);}
if (!ReadXMLData(J,xsm::Tz,xSystem)) {returnstring<<"T"<<int(J);}
if (!ReadXMLData(readstring,xsm::P,xSystem)) {returnstring<<readstring;}
return returnstring.str();
}
}
#endif //__XLABEL__CXX__

