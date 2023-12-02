/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file davidson_vfile.cpp
\brief Davidson diagonalization of matrix in the file. Matrix is loaded in the memory for processing
*/
#include <iostream>
#include <av/argparser.cxx>
#include <av/statusclass.cxx>
#define CSTATUS
#include <xml/xmlParser.cxx> //needed only if input is xml
using xml::XMLNode;
#include <xsm/schema.h>
#include <cstdlib>
#include <tav/sparse.cxx>
#include <kdav/davidson.cxx>
using std::atoi;




int main (int argc, char* argv[]) {
av::argparser args(argc,argv);
args.SetDescription("Davidson diagonalization of a matrix stored in a file; matrix is loaded in the memory for processing.");
args.SetDefault("File where matrix to be diagonalized is stored file.HH",2,2);
args.SetFlag("i","iterations","Set the maximum number of iterations.",0,2);
args.SetFlag("n","number","Set the desired number of states.",0,2);
args.SetFlag("b","block","Set the block size.",0,2);
args.SetFlag("f","file","File name where initial set of fectors is stored",0,2);
if(args.ParseArguments()) return 1;
string sysfile=args(1);

 //all below is done in case file is a matrix in xml
// if (sysfile.back() != 'H') {
 if (sysfile[sysfile.length()-1] != 'H') {
    XMLNode xMainNode=XMLNode::openFileHelper(args(1).c_str(), xsm::cosmo);
    if (xMainNode.isEmpty()) FatalError("No data in xml");
    XMLNode xSystemNode=xml::SelectNode(xMainNode,xsm::system,"name");
    XMLNode xMatrixNode=xSystemNode.getChildNode(xsm::system_matrix,-1); //get last child
    if (xMatrixNode.isEmpty()) FatalError("You need Hamiltonian matrix");
   xml::AskAttribute(sysfile,"name", "Matrix file name: ", xMatrixNode);
    cerr<<"Diagonalizing file "<<sysfile<<endl; 
 }


int iters=0; int bloc=10; int nstate=0;
    string invectorf; //file for in vectors
if (args.isSet("i",1)) iters=args.intarg("i",1);
if (args.isSet("n",1)) nstate=args.intarg("n",1);
if (args.isSet("b",1)) bloc=args.intarg("b",1);
if (args.isSet("f",1)) { invectorf=args("f",1); }
//cerr<<iters<<" "<<nstate<<" "<<bloc<<endl;
if (bloc<nstate) bloc=nstate; //block must be bigger then number of saved states

cerr<<"Work with: "<<sysfile<<endl;

tav::SparseMatrix<ftyp> FHH;
FHH.Read(sysfile.c_str());
FHH.type=1;
//remember this thing needs a pointer 

//output file
string outf=sysfile+".EE";
//DavidsonDiagonalization(outf, FHH, bloc,iters,nstate);
if (args.isSet("f",1))  DavidsonDiagonalization(outf.c_str(), FHH,invectorf.c_str(), bloc,iters,nstate);
   else DavidsonDiagonalization(outf.c_str(), FHH, bloc,iters,nstate);
return 0;
		
}

