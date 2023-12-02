/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


#include <debug.h> //AskUser
#include <xml/xmlParser.cxx>
#include <sm/spinclass.cxx> //spin
#include <sm/spherical2sp.cxx> //SPTPlus
#include <sm/EMspherical2sp.cxx>
#include <rw/basic-rw.cxx>
//using rw::TestFile


XMLNode MakeXMLOperator(XMLNode xMainNode) {
    
    char ntype;
    AskUser("Operator type Electric/Magnetic/GamowTeller/Fermi (E/M/G/F) ",ntype);
    int K;
    AskUser("Enter Multipolarity ",K);
    XMLNode xOperator;
    if (ntype=='E') {
        double en,ep; //effective charges
        cerr<<"We are using electric transition E"<<K<<" \n";
        if (K==1) {
            cerr<<"This is a dipole transition "<<K<<" \n";
            try {
                int A;
                XMLNode xsystemnode=xMainNode.getChildNode(xsm::system);
//                getXMLData(A,xsm::nucleus_A,xMainNode);
                getXMLData(A,xsm::nucleus_A,xsystemnode);
                SM::spin Tz;
//                getXMLData(Tz, xsm::Tz, xMainNode);
                getXMLData(Tz,xsm::Tz,xsystemnode);
                int NZ=int(Tz);
                int N=(NZ+A)/2;
                int Z=(A-NZ)/2;
                en=-double(Z)/double(A);//From Ring Schuck Appendix B p 592
                ep=double(N)/double(A);//From Ring Schuck Appendix B p 592
                cout<<"Neutron effective charge: "<<en<<endl;
                cout<<"Proton effective charge: "<<ep<<endl;

            }
            catch (...)
            {
                AskUser("Enter neutron effective charge: ", en);
                AskUser("Enter proton effective charge: ", ep);
            }
        }
        else {
            AskUser("Enter neutron effective charge: ", en);
            AskUser("Enter proton effective charge: ", ep);
        }
        xOperator=xMainNode.addChild(xsm::system_operator, false, 1);
        char qname[20]; sprintf(qname,"%c%i",ntype,K);
        // cerr<<qname<<endl;
        addAttribute(xOperator,"name", qname);
        addAttribute(xOperator,xsm::J,K);
        addAttribute(xOperator,"en",en);
        addAttribute(xOperator,"ep",ep);
        //radial part
        char rfile[256];
        AskUser("Enter file for radial overlap (press return to use HO) ",rfile);
        if (!rw::TestFile(rfile)) {
            cerr<<"We will use HO with hw=1 \n";
            /*
             cerr<<"Using HO \n"
             int A;
             try {getXMLData(A,xsm::nucleus_A,xSystem);}
             catch (...) {AskUser("Enter mass number to use for radial scaling A", A);}
             */
        } //end file does not exist
        else {
            XMLNode xRadial=xOperator.addChild("radial");
            addAttribute(xRadial,"file",rfile);
            addAttribute(xRadial,"scale",K);
        }
    } //end electric
    if (ntype=='M') {
        double mn,mp; //magnetic moments
        //	AskUser("Enter neutron effective charge: ", mn);
        //	AskUser("Enter neutron effective charge: ", mp);
        xOperator=xMainNode.addChild(xsm::system_operator, false, 1);
        char qname[20]; sprintf(qname,"%c%i",ntype,K);
        // cerr<<qname<<endl;
        addAttribute(xOperator,"name", qname);
        addAttribute(xOperator,xsm::J,K);
        //addAttribute(xOperator,"mn",en);
        //addAttribute(xOperator,"mp",ep);
    }
    if (ntype=='G') {
        xOperator=xMainNode.addChild(xsm::system_operator, false, 1);
        // cerr<<qname<<endl;
        int tflag; //tflag=1 increase isospin beta+
        AskUser("Beta + or beta - operator (enter +1 or -1) ", tflag);
        if (tflag==1) addAttribute(xOperator,"name", "GT+");
        else addAttribute(xOperator,"name", "GT-");
        //addAttribute(xOperator,"type", tflag);
        addAttribute(xOperator,xsm::J,K);
        //addAttribute(xOperator,"mn",en);
        //addAttribute(xOperator,"mp",ep);
    }
    if (ntype=='F') {
        xOperator=xMainNode.addChild(xsm::system_operator, false, 1);
        // cerr<<qname<<endl;
        addAttribute(xOperator,"name", "F");
        addAttribute(xOperator,xsm::J,0);
        //addAttribute(xOperator,"mn",en);
        //addAttribute(xOperator,"mp",ep);
    }
    return xOperator;
}

int ReadXMLOperator(
                    MBOperator &MUME, ///<[o] Output two-body operator
                    int kappa, ///<[i] Input projection
                    SM::SingleParticleStates &sps, ///< [i] single-particle states
                    XMLNode xOperator, ///< [i] xml pointer to operator
                    bool phaseConvention = true
)
{
    if (xOperator.isEmpty()) CriticalError("Operator is not available in xml");
    //count from sps how many levels we have
    int levels=0;
    for (int s=0;s<sps.size();s++) if (levels<sps[s].level) levels=sps[s].level;
    levels++;
    //
    tav::Tensor<2,double> RE(levels,levels); //even multipole
    ifstream rfile;
    //for (int i=0;i<x.size();i++) for (int j=0;j<x.size();j++)  RE[i][j]=1.0;
    XMLNode xRadial=xOperator.getChildNode("radial");
    if (!xRadial.isEmpty()) {
        string rfilename;
        getXMLData(rfilename,"file",xRadial);
        rfile.open(rfilename.c_str());
        rfile>>RE;
        rfile.close();
        cout<<RE<<endl;
        Pause("Verify radial overlaps, enter to continue");
    }
    
    int K=2;
    getXMLData(K,xsm::J,xOperator);
    char tname[256];
    getXMLData(tname,"name",xOperator);
    //Pause(tname);
    //MBOperator MUME(2);
    if (tname[0]=='G') //beta
    {
        
        //if (PT.Tz = DT.Tz) FatalError("Isospin must change in beta decay");
        //Must be beta transition
        int ttype;
        //	AskAttribute(ttype,"type","We need beta operator type (+1 or -1)? :", xOperator);
        if (tname[2]=='+') ttype=1; else ttype=-1;
        if (K!=1) FatalError("We only use vector transitions for beta decay");
        if (!xRadial.isEmpty()) BetaJAMultipole(MUME,sps, K, kappa, RE, 0);
        else BetaJAMultipole(MUME,sps, K, kappa, 0, ttype);
        //   cerr<<MUME;
        //   Pause();
        MUME*=sqrt(4.0*PI); //change to cartesian so the units are g^2/4\pi
    } //end beta decay
    
    if (tname[0]=='F') //beta
    {
        
        //if (PT.Tz = DT.Tz) FatalError("Isospin must change in beta decay");
        //Must be beta transition
        if (K!=0) FatalError("Fermi operator is scalar");
        SM::SPTPlus(MUME, sps);
        SM::Hermitian2FullOperator(MUME);
        //   cerr<<MUME;
        //   Pause();
    } //end beta decay
    
    if (tname[0]=='E') //electric
    {
        //if ((K&1)^(PT.P)^(DT.P)) FatalError(xsm::E<<K<<" oprator does not change parity");
        double en;
        AskAttribute(en,"en","Neutron effective charge :", xOperator);
        double ep;
        AskAttribute(ep,"ep","Proton effective charge :", xOperator);
        if (!xRadial.isEmpty()) ElectricSPMultipole(MUME, sps, K, kappa, RE, en, ep);
        else ElectricSPMultipole(MUME, sps, K, kappa, en, ep,phaseConvention);
    }
    
    if (tname[0]=='M') //magnetic
    {
        double mun=phy::muneutron;
        InquireAttribute(mun,"mu_n","Neutron effective magnetic moment :", xOperator);
        double mup=phy::muproton;
        InquireAttribute(mup,"mu_p","Proton effective magnetic moment :", xOperator);
        double gl_n=0.0;
        InquireAttribute(gl_n,"gl_n","Neutron orbital gyromagnetic ratio :", xOperator);
        double gl_p=1.0;
        InquireAttribute(gl_p,"gl_p","Proton orbital gyromagnetic ratio :", xOperator);
        
        
        
        if (!xRadial.isEmpty()) MagneticSPMultipole(MUME, sps,  K, kappa,RE,mun,mup,gl_n,gl_p);
        MagneticSPMultipole(MUME, sps,  K, kappa,mun,mup,gl_n,gl_p,phaseConvention);
    }   
    
    //operator is done
    return K;
} 		

XMLNode InquireXMLOperator(XMLNode xMainNode) {
    XMLNode xOperator=xMainNode.getChildNode(xsm::system_operator);
    if (xOperator.isEmpty()) {
        MakeXMLOperator(xMainNode);
        xOperator=xMainNode.getChildNode(xsm::system_operator);
    }
    else { 
        if (Inquire("Use avaliable operators ",true)) xOperator=xml::SelectNode(xMainNode,xsm::system_operator,"name");
        else {xOperator=MakeXMLOperator(xMainNode); }
    }
    return xOperator;
}
