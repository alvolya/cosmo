/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file DXHV.cxx
\brief Davidson diagonalization of Virtual Hamiltonian (generated on fly) 
\ingroup gp_xsm
*/
#ifndef __DXHV__CXX__
#define __DXHV__CXX__

#include <Eigen/Eigenvalues>
using namespace Eigen;

#include <xsm/xManyBodyStates.cxx>
#include <xsm/xsphericalint.cxx>
#include <sm/MakeSPO.cxx>
#include <sm/spherical2sp.cxx>
#include <tav/sparse.cxx>
#include <tav/csbmatrix.cxx>
#include <xsm/xReadSPInt.cxx>
#include <sm/PPMatrix.cxx>
#include <kdav/davidson.cxx>

namespace xsm{

    int DXHV(XMLNode xSystemNode, XMLNode xMainNode, int iters, int bloc, int nstate, string invectorf,int limit = 1000 ) {
 SM::SingleParticleStates sps; Make(sps,xMainNode.getChildNode(xsm::valencespace)); 
 string sysname;
 xml::AskAttribute(sysname,"name", "Enter system name: ", xSystemNode);
 //many-body states
 SM::Many_Body_States st;
SM::Get(st,xMainNode,xSystemNode);
//hamiltonian node
string hname=sysname+".HH";
//hamiltonian node
tav::STensor<spsint,double> ESP(2); //output one-body interaction
 tav::STensor<spsint,double> VSPX(4);
 
 xReadSPInt(ESP, VSPX, xSystemNode, xMainNode,sps);
   
  SM::PPMatrix  HH(st);
  HH.push_back(&ESP);
  HH.push_back(&VSPX);
  //cerr<<VSPX<<endl;
   //XMLNode xNode=xMainNode.addChild("Hamiltonian");
   //SM::WriteXMLPPMatrix(xNode,HH); 
  // cout<<xMainNode.createXMLString(true)<<endl;
 //Pause(); 


  
  cerr<<"generating "<<hname<<endl;
  cerr<<"dimension "<<st.n<<endl;
    if (st.n>limit && limit>0)
    {
//        tav::SparseMatrix<ftyp> FHH(st.n);
//        PPMatrix2SparseMatrix(FHH,HH);
        uint16_t br = std::sqrt(1.*st.n) + 1;
        tav::SparseMatrixCSB<double,uint32_t,uint16_t> FHH(st.n,br,br,true);
        PPMatrix2CSBBlk(FHH,HH);
        hname+=".EE";
 // DavidsonDiagonalization(hname.c_str(), FHH, bloc,iters,nstate);
        if (rw::TestFile(invectorf.c_str())) DavidsonDiagonalization(hname.c_str(), FHH,invectorf.c_str(), bloc,iters,nstate);
        else DavidsonDiagonalization(hname.c_str(), FHH, bloc,iters,nstate);
    }
    else
    {
        std::cout << "Matrix is small for davidson. Using exact diagonalization" << std::endl;
        MatrixXd FHH=MatrixXd::Zero(st.n,st.n);
        PPMatrix2MatrixXd(FHH,HH);
        hname+=".EE";
        
        //-----------THIS IS COPY FROM EXACTEV------------------------------------------------------------
        Eigen::SelfAdjointEigenSolver<MatrixXd> eis(FHH);
        
        for (int i=0;i<nstate;i++)
            cout<<eis.eigenvalues()[i]<<" ";
        if (nstate==0)
        {
            for (int i=0;i<st.n;i++)
                cout<<eis.eigenvalues()[i]<<" ";
        }
        cout<<endl;
        
        ofstream outf;
        outf.open(hname.c_str(), std::ios::binary);
        outf.write((char*)(&(st.n)),sizeof(unsigned)); //save number of eigenvalues
        // cout<<EE<<endl;
        double zerror=0.0; //write error as zero
        for (int ss=0;ss<st.n;ss++)
        {
            outf.write((char*)(&(eis.eigenvalues()[ss])),sizeof(double)); //save eigenvalues
            outf.write((char*)(&zerror),sizeof(double)); //save error in eigenvalue
            outf.write((char*)(&zerror),sizeof(double)); //save error in eigenvector
            outf.write((char*)(&st.n),sizeof(uint_nbasis)); //save dimension
            VectorXd vv=eis.eigenvectors().col(ss);
            //     double *prt=new double[n];
            //     prt=&(eis.eigenvectors().col(ss)[0]);
            //     delete[] prt;
            outf.write((char*)(&(vv[0])),st.n*sizeof(double)); //save eigenvector
        }
        //-----------END COPY FROM EXACTEV---------------------------------------------------------------

    }
  //HH.Write(hname.c_str());
return 0;
}

} //end namespace 






#endif //__DXHV__CXX__

