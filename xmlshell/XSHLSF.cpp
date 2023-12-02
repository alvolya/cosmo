/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/



/*!
 \author Alexander Volya  <http://www.volya.net>
 \file XSHLSF.cpp
 \brief Spectroscopic factors and amplitudes
 
 Usage: XSHLSF parent.xml daughter.xml
 if number of particles is the same we do an overlap
 for number of particles different by one we do SF
 */

#if __cplusplus <= 199711L
typedef unsigned short spsint; //this part is for old c++98
#else
#include <cstdint>
typedef uint16_t spsint;
#endif
typedef unsigned long long  uint_nbasis;
    //#include <math.h>
#define INFO yes
    //#include <debug.h>
    //#include <tav/sparse.cxx>
    //#include <av/ThreeJSymbol.h>
    //#include <sm/ModelSpace.cxx>
#include <av/argparser.cxx>
    //#include <sm/SingleParticleStates.cxx>
#include <xsm/xManyBodyStates.cxx>
#include <sm/MakeSPO.cxx>
#include <xsm/xState.cxx>
#include <sm/EMspherical2sp.cxx>
#include <tav/ThreeJSymbolM.cxx>
using tav::ThreeJSymbol;
#include <sm/PPMatrix.cxx>
    //#include <sm/supplemental.cxx>


/*-----------MAIN---------------------------------------*/
int main(int argc, char** argv) {
    av::argparser args(argc,argv);
    args.SetDescription("Analyze many-body transition operator");
    args.SetDefault("XML file containing a parent system description.",2,2);
    args.SetFlag("f","final","XML file containg a final daughter system, use if final system is in a separate xml",0,2);
    args.SetFlag("o","output","Output file with plain text",0,2);
    args.SetFlag("of","output-full","Output file with plain text and full data that includes j",0,1);
    args.SetFlag("s","square","Compute operator's square (non-hermitian)",0,1);
    if(args.ParseArguments()) return 1;
    
    string sysfile=args(1);
    XMLNode xMainNode=XMLNode::openFileHelper(sysfile.c_str(), xsm::cosmo); //we start by opening xml input file
    
    
    
    XMLNode xValence=xMainNode.getChildNode(xsm::valencespace);
    int levels=xValence.nChildNode(xsm::valencespace_orbital);
        //cerr<<levels<<endl;
    SM::SingleParticleStates sps; Make(sps,xValence);
    XMLNode xParentSys=xml::SelectNode(xMainNode,xsm::system,"name");
        //XMLNode xParentSys=xMainNode.getChildNode(xsm::system,-1); //get last child
        //if (xParentSys.isEmpty()) {FatalError("system must exit"); }
    char sysname[256];
    xml::AskAttribute(sysname,"name", "Matrix file name: ", xParentSys);
    cerr<<"Parent system "<<sysname<<endl;
        //ReadXMLQN(x,xSystemNode);
        //many-body states
    SM::Many_Body_States stp;
    SM::Get(stp,xMainNode,xParentSys);
    
    
        //daughter xml
    XMLNode xMainNodeDaughter=xMainNode; // daughter system
    if (args.isSet("f",1)) {
        xMainNodeDaughter=XMLNode::openFileHelper(args("f",1).c_str(), xsm::cosmo);
    }
    cerr<<"Select daughter system \n";
    XMLNode xDaugterSys=xml::SelectNode(xMainNodeDaughter,xsm::system,"name");
    
    xml::AskAttribute(sysname,"name", "Matrix file name: ", xDaugterSys);
    cerr<<"Daughter system "<<sysname<<endl;
    SM::Many_Body_States std;
    SM::Get(std,xMainNodeDaughter,xDaugterSys);
        //SM::Make(std,xMainNoded);
    
    
    
    
    
    SM::State PT, DT;
    ReadXMLQN(PT, xParentSys);
    getXMLData(PT.N ,xsm::valence_A, xParentSys);
    ReadXMLQN(DT, xDaugterSys);
    getXMLData(DT.N ,xsm::valence_A, xDaugterSys);
    
    
    
    int pstates=xParentSys.nChildNode(xsm::system_eigenstate);
    if (pstates==0) FatalError("no states in xml");
    int dstates=xDaugterSys.nChildNode(xsm::system_eigenstate);
    if (dstates==0) FatalError("no states in xml");
    string stmp;
    XMLNode xParentStateNode;
    XMLNode xDaugterStateNode;
        // int Jz=PT[0].Jz-DT[0].Jz;
    std::ofstream ofile;
    if (args.isSet("o",1)) {ofile.open(args("o",1).c_str());
        ofile<<"#Spectroscopic factors for transition between ";
        ofile<<getXMLData(stmp,xsm::name,xParentSys)<< " and ";
        ofile<<getXMLData(stmp,xsm::name,xDaugterSys)<<endl;
        ofile<<"#File: "<<sysfile<<" ";
        if (args.isSet("f",1)) ofile<<"Final(daugter) system file: "<<args("f",1);
        ofile<<endl;
        if (!args.isSet("of",1)) ofile<<"#parent Ex(p) daugther Ex(d) L SF=|A|^2 \n";
            else ofile<<"#parent Ex(p) daugther Ex(d) L J A SF=|A|^2 \n";
    }
    
    
    double tsum=0.0;
    
    for (int p = 0; p < pstates; p++) {
        xParentStateNode=xParentSys.getChildNode(xsm::system_eigenstate,p);
        ReadXMLState(PT, xParentStateNode); //read the rest of qn
        string ParentName;
        getXMLData(ParentName,"name",xParentStateNode);
        
        for (int q = 0; q < dstates; q++) {
            xDaugterStateNode=xDaugterSys.getChildNode(xsm::system_eigenstate,q);
                //	cout<<xDaugterStateNode.createXMLString(true);
            ReadXMLState(DT, xDaugterStateNode); //read the rest of qn
            string DaugterName;
            getXMLData(DaugterName,"name",xDaugterStateNode);
            
            
            
            int Jz = PT.Jz - DT.Jz;
            int Tz = PT.Tz - DT.Tz;
            bool P=PT.P^DT.P; //parity
                              //cout << PT << DT;
            cout<<"Parent: "<<ParentName<<" Ex="<<getXMLData(stmp,xsm::Ex,xParentStateNode)<<" E="<<PT.E<<endl;
            cout<<"Daughter: "<<DaugterName<<" Ex="<<getXMLData(stmp,xsm::Ex,xDaugterStateNode)<<" E="<<DT.E<<endl;
                // cerr<<PT.N<<" "<<DT.N<<endl;
            if ((PT.N - DT.N) == 0) {
                    //doing overlap
                if (Jz!=0) continue;
                if (Tz!=0) continue;
                if (P) continue;
                if (PT.J!=DT.J) continue;
                if (PT.T!=DT.T) continue;
                double tmp_overlap=DT.z*PT.z;
                cout<<"<x|y> "<<tmp_overlap<<" |<x|y>|^2 = "<<tmp_overlap*tmp_overlap<<" "<<endl;
                continue;
            }
            
            if ((PT.N - DT.N) != 1) FatalError("The number of particles in final system must be less by one than in initial");
            int Lmax=20;
            double L_SF[Lmax]; //array for L-values
            for (int i=0;i<Lmax;++i) L_SF[i]=0; //initialize array of max
            for (spsint s = 0; s < sps.size(); s++) {
                if (sps[s].Jz != Jz) continue;
                if (sps[s].Tz != Tz) continue;
                if (sps[s].P^P) continue;
                double htmp = tav::ClebschGordan(DT.J, DT.Jz, sps[s].J, Jz, PT.J, PT.Jz);
                if (tav::Abs(htmp) < 1E-10) continue;
                    // At this stage we can compute transition
                /*! \note It has been tested that the answer does not change if complex conjugated, annihilation operator is evaluated.
                 */
                    // the following lines also work fine
                    //        #include <sm/Den1B.cxx>
                    //         tav::Tensor<1,double> X(stp.n);
                    //        SM::TransitionAmplitude(s, DT.z.elem, std, X.elem,  stp) ;
                    //        double ovtmp1 = X*PT.z;
                SM::PPMatrix A(std, stp);
                MBOperator a(1);
                a(s) = 1;
                    //	cout<<"s= "<<int(s)<<" "<<sps[s].J<<" "<<" Sjz= "<<sps[s].Jz<<" Stz="<<sps[s].Tz<<" Jz= "<<Jz<<" Tz= "<<Tz<<endl;
                A.push_back(&a);
                double ovtmp = A.OverlapNH(DT.z, PT.z);
                cout << sps[s].tag << "/2 C2S=" << tav::Sqr(ovtmp / htmp) << " A=" << (ovtmp / htmp) << endl;
                           if (args.isSet("of")&&args.isSet("o",1)) {
                             ofile<<getXMLData(stmp,xsm::name,xParentStateNode)<<"\t";
                             ofile<<getXMLData(stmp,xsm::Ex,xParentStateNode)<<"\t";
                             ofile<<getXMLData(stmp,xsm::name,xDaugterStateNode)<<"\t";
                             ofile<<getXMLData(stmp,xsm::Ex,xDaugterStateNode)<<"\t";
                             ofile<<sps[s].L<<"\t";
                             ofile<<sps[s].J<<"\t";
                             tsum+=tav::Sqr(ovtmp / htmp);
                             ofile<< (ovtmp / htmp) << "\t "<<tav::Sqr(ovtmp / htmp) << endl;
                             }
                L_SF[sps[s].L]+=tav::Sqr(ovtmp / htmp);
                
            } //s.p. sf loop
            cout<<endl;
            if (args.isSet("o",1)&&(!args.isSet("of")))
            for (int i=0;i<Lmax;++i) if (L_SF[i]>0) {
                ofile<<getXMLData(stmp,xsm::name,xParentStateNode)<<"\t";
                             ofile<<getXMLData(stmp,xsm::Ex,xParentStateNode)<<"\t";
                             ofile<<getXMLData(stmp,xsm::name,xDaugterStateNode)<<"\t";
                             ofile<<getXMLData(stmp,xsm::Ex,xDaugterStateNode)<<"\t";
                             ofile<<i<<"\t";
                    //ofile<<sps[s].J<<"\t";
                    //       tsum+=tav::Sqr(ovtmp / htmp);
                             ofile<< L_SF[i]  << endl;
            } //end loop over L
            
        } //daughter loop
    } //parent loop
    return 0;
}









