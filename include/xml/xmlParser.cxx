/*! \file xmlParser.cxx
\ingroup gp_nr
*/
#ifndef __XMLPARSER__CXX__
#define __XMLPARSER__CXX__
#include <iostream>
#include <sstream>
#include <debug.h>
#include <cstdlib> //atoi atof functions
#include <cstdio> //sprintf
#include <inquire.cxx>
using tav::Inquire;
namespace xml{
#include "xmlParser.h"
#include "xmlParser.cpp"

/*! \brief Select a childnode with a given selector name
*/
XMLNode SelectNode(XMLNode xMainNode, const char* nodename, const char* atr_selector) {
int nsp=xMainNode.nChildNode(nodename);
if (nsp==0) CriticalError("Node "<<nodename<<" is not found\n");
if (nsp==1) return xMainNode.getChildNode(nodename);
XMLNode xNodeSelect;
std::vector<std::string> namelist;
std::cerr<<"Avaialble "<< nodename<<"'s :\n";
for (int i=0;i<nsp;i++) {
xNodeSelect=xMainNode.getChildNode(nodename,i);
namelist.push_back(xNodeSelect.getAttribute(atr_selector)); //get this node in the list
int it=xNodeSelect.nAttribute(); //get number of attributes
std::cerr<<atr_selector<<" "<<xNodeSelect.getAttribute(atr_selector)<<std::endl;
for (int atr=0;atr<it;atr++) 
std::cerr<<xNodeSelect.getAttributeName(atr)<<"="<<xNodeSelect.getAttributeValue(atr)<<"  ";
std::cerr<<std::endl;
std::cerr<<std::endl;
}
char  myselection[256];
XMLNode xReturn;
do {
/*
std::cerr<<"Select "<< atr_selector<<" :";
std::cin>>myselection;
//AskUser("select: ", myselection);
*/
xReturn=xMainNode.getChildNodeWithAttribute(nodename,atr_selector,Inquire(namelist).c_str());
if (xReturn.isEmpty()) std::cerr<<"The node with "<<atr_selector<<"="<<myselection<<" not found, try again \n";
} while (xReturn.isEmpty());
return xReturn;
}


/*! \brief Trim main node a childnode with a given selector name
*/
int SelectNodes(
XMLNode xMainNode, ///<[i/o] node to change
const char* nodename, ///<[i] name of child node
bool bdef=false ///<default selection
) {
int nsp=xMainNode.nChildNode(nodename);
XMLNode xNodeSelect;
if (nsp==0) CriticalError("Node "<<nodename<<" is not found\n");
int i=0;
std::cerr<<std::endl;
while (i<nsp) {
xNodeSelect=xMainNode.getChildNode(nodename,i);
int it=xNodeSelect.nAttribute(); //get number of attributes
for (int atr=0;atr<it;atr++) 
std::cerr<<xNodeSelect.getAttributeName(atr)<<"="<<xNodeSelect.getAttributeValue(atr)<<"  ";
std::cerr<<std::endl;
if (!Inquire("Select this node ?",bdef)) {xNodeSelect.deleteNodeContent(); nsp--;}
else i++;
}; 
return 1;
}


//add via string stream
template <class cnumber> 
void updateAttribute( 
XMLNode xNode, ///<node 
const char* xmlname, ///< name of xml parameter
const cnumber &x
) {
if (xNode.isEmpty()) CriticalError("xml node is empty");
//try to find data as attribute
std::stringstream ss (std::stringstream::in | std::stringstream::out);
const char *t=xNode.getAttribute(xmlname); //try to get parameter
if (t!=NULL) xNode.deleteAttribute(xmlname); //attribute is present
ss<<x;
char tmp[256];
ss>>tmp;
xNode.addAttribute(xmlname,tmp);
}

//copy all atrigutes
void copyAttributes(XMLNode outxNode, XMLNode xNode)
{
    if (xNode.isEmpty()) CriticalError("xml input node is empty");
    if (outxNode.isEmpty()) CriticalError("xml output node is empty");
    int it=xNode.nAttribute(); //get number of attributes
    for (int atr=0;atr<it;atr++){
        updateAttribute(outxNode,xNode.getAttributeName(atr),xNode.getAttributeValue(atr));
//std::cerr<<xNode.getAttributeName(atr)<<"="<<xNode.getAttributeValue(atr)<<"  ";
    }//end loop over attributes
} //end void function


//add attribute via string stream
template <class cnumber> 
void addAttribute( 
XMLNode xNode, ///<node 
const char* xmlname, ///< name of xml parameter
const cnumber &x
) {
if (xNode.isEmpty()) CriticalError("xml node is empty");
//try to find data as attribute
std::stringstream ss (std::stringstream::in | std::stringstream::out);
const char *t=xNode.getAttribute(xmlname); //try to get parameter
if (t!=NULL) CriticalError("addAttribute: Attribute "<<xmlname<<"is already present"); //attribute is present
ss<<x;
char tmp[256];
ss>>tmp;
xNode.addAttribute(xmlname,tmp);
}


/** \breif Using stringstream get arbitrary xml data from Attribute or Child
*/
template <class cnumber> 
cnumber& getXMLData(cnumber &x, 
const char* xmlname, ///< name of xml parameter
XMLNode xNode ///<node 
) {
if (xNode.isEmpty()) CriticalError("xml node is empty");
//try to find data as attribute
std::stringstream ss (std::stringstream::in | std::stringstream::out);
const char *t=xNode.getAttribute(xmlname); //try to get parameter
if (t!=NULL) {
ss<<t;
ss>>x;
return x;
}
//try to find data as child
XMLNode xTNode=xNode.getChildNode(xmlname);
if (!xTNode.isEmpty()) {
ss<<xTNode.getText();
ss>>x;
return x;
}
CriticalError("data "<<xmlname<<" is not found ");
}


/** \breif Using stringstream read arbitrary xml data return 1 if fail
*/
template <class cnumber> 
int ReadXMLData(cnumber &x, 
const char* xmlname, ///< name of xml parameter
XMLNode xNode ///<node 
) {
if (xNode.isEmpty()) return 1; //CriticalError("xml node is empty");
//try to find data as attribute
std::stringstream ss (std::stringstream::in | std::stringstream::out);
const char *t=xNode.getAttribute(xmlname); //try to get parameter
if (t!=NULL) {
ss<<t;
ss>>x;
return 0;
}
//try to find data as child
XMLNode xTNode=xNode.getChildNode(xmlname);
if (!xTNode.isEmpty()) {
ss<<xTNode.getText();
ss>>x;
return 0;
}
return 1; // CriticalError("data "<<xmlname<<" is not found ");
}


/*
This reads xml data but returns the default if error or if not found
*/
template <class cnumber> 
cnumber& InquireXMLData(cnumber &x, //default 
const char* xmlname, ///< name of xml parameter
XMLNode xNode ///<node 
) {
if (xNode.isEmpty()) return x;
//try to find data as attribute
std::stringstream ss (std::stringstream::in | std::stringstream::out);
const char *t=xNode.getAttribute(xmlname); //try to get parameter
if (t!=NULL) {
ss<<t;
ss>>x;
return x;
}
//try to find data as child
XMLNode xTNode=xNode.getChildNode(xmlname);
if (!xTNode.isEmpty()) {
ss<<xTNode.getText();
ss>>x;
return x;
}
//CriticalError("data "<<xmlname<<" is not found ");
return x;
}


//here we provide some additional xml function 
//use string stream to Ask= means check xml and then ask user and update
/*! \breif This function will read from xml or ask user and add to xml text data in node
<somedata>12.5</somedata>
*/
template <class cnumber> 
cnumber& AskNode(
cnumber &x, ///< variable 
const char* xmlname, ///< name of xml parameter
const char* question, ///< question to ask from user
XMLNode xNode ///<node 
) {
std::stringstream ss (std::stringstream::in | std::stringstream::out);
XMLNode xTNode=xNode.getChildNode(xmlname);
//add text to node
if (xTNode.isEmpty()) {
xTNode=xNode.addChild(xmlname);
//AskUser(x,question);
std::cerr<<question;
std::cin>>x;
ss<<x;
char text[50];
ss>>text;
xTNode.addText(text);
return x;
}
ss<<xTNode.getText();
ss>>x;
return x;
}


/*! \breif This function will read from xml or ask user and add to xml text data in node
<node somedata=12.5>
*/
template <class cnumber> 
cnumber& AskAttribute(
cnumber &x, ///< variable 
const char* xmlname, ///< name of xml parameter
const char* question, ///< question to ask from user
XMLNode xNode ///<node 
) {
std::stringstream ss (std::stringstream::in | std::stringstream::out);
const char *t=xNode.getAttribute(xmlname); //try to get parameter 
//add text to node
if (t==NULL) {
std::cerr<<question;
std::cin>>x;
//cerr<<" read x "<< x<<endl;
ss<<x;
char text[50];
ss>>text;
//cerr<<" read x(text) "<< text<<endl;
xNode.addAttribute(xmlname,text);
return x;
}
ss<<t;
ss>>x;
return x;
}

/*! \breif This function will read from xml or ask user and add to xml text data in node
<node somedata=12.5>
*/
template <class cnumber> 
cnumber& InquireAttribute(
cnumber &x, ///< variable 
const char* xmlname, ///< name of xml parameter
const char* question, ///< question to ask from user
XMLNode xNode ///<node 
) {
std::stringstream ss (std::stringstream::in | std::stringstream::out);
const char *t=xNode.getAttribute(xmlname); //try to get parameter 
//add text to node
if (t==NULL) {
Inquire(x,question);
//cerr<<" read x "<< x<<endl;
ss<<x;
char text[50];
ss>>text;
//cerr<<" read x(text) "<< text<<endl;
xNode.addAttribute(xmlname,text);
return x;
}
ss<<t;
ss>>x;
return x;
}

/*!
\brief Delete all xml subnodes with this name
*/
int DeleteXMLNodes(
XMLNode &xOutMainNode, ///<main
const char* nname ///<name to delete
)
{
 int nsystem=xOutMainNode.nChildNode(nname);
        if (nsystem!=0) {
        std::cerr<<"The nodes \""<<nname<<"\" are about to be deleted, ";
        if (Inquire("delete old nodes in xml ",false)) {
        for (int in=0;in<nsystem;in++) (xOutMainNode.getChildNode(nname)).deleteNodeContent();
}
}
return 0;
}

}

#endif //__XMLPARSER__CXX__
