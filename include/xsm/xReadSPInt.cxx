/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! \file xReadSPInt.cxx
 \ingroup gp_xsm
 \author Alexander Volya  <http://www.volya.net>
 \file XHH+JJ.cpp
 \brief Full reader of xml interactions xml input, sp-interaction output
 */
#ifndef __XREADSPINT__CXX__
#define __XREADSPINT__CXX__
#include <xsm/xsphericalint.cxx>
#include <sm/MakeSPO.cxx>
#include <sm/spherical2sp.cxx>
#include <sm/ReadSPInt.cxx>
#include <kelem/CM_Lawson.cxx>
#include <kelem/rk_Interaction.cxx>
#include <kelem/delta.cxx>
#include <sm/AddQQVspJT.cxx>
int xReadSPInt(
               tav::STensor<spsint,double> &ESP, //output one-body interaction
               tav::STensor<spsint,double> &VSPX, //output two-body interaction
               XMLNode xSystemNode, ///< [i/o] system node
               XMLNode xMainNode, ///< [in] need for everything
               SM::SingleParticleStates &sps ///< [in] single particle states
)
{
    
    //figure out J^2 term
    string sysname;
    xml::AskAttribute(sysname,"name", "Enter system name: ", xSystemNode);
    string hname=sysname+".HH";
    
    
    XMLNode xMatrixNode=xSystemNode.getChildNode(xsm::system_matrix);
    if (xMatrixNode.isEmpty()) {
        cerr<<"Creating matrix node"<<endl;
        xMatrixNode=xSystemNode.addChild(xsm::system_matrix);
        addAttribute(xMatrixNode,"name",hname);
        // jadd=0.0;
        //  tav::Inquire(jadd,"Enter scaling for J+J- part ");
        //  addAttribute(xMatrixNode,"JJ",jadd);
        addAttribute(xMatrixNode,"format","pp");
    }
    else {
        updateAttribute(xMatrixNode,"name",hname);
        updateAttribute(xMatrixNode,"format","pp");
    }
    
    
    
    int N; //number of particles
    AskAttribute(N,xsm::valence_A,"Enter valence particle number: ",xMainNode.getChildNode(xsm::system));
    XMLNode xHamiltonianNode=xMainNode.getChildNode(xsm::hamiltonian) ;
    cerr<<"reading interaction files "<<endl;
    if (xHamiltonianNode.isEmpty()) FatalError("ReadSMinINT:  empty node");
    int files=xHamiltonianNode.nChildNode("file");
    int tagfiles = xHamiltonianNode.nChildNode("tagfile");
    
    if (files!=0) { //if there are files then this is oxbash format
        cerr<<"Interaction files are in OXBASH format "<<endl;
        //_____________________INTERACTIONS__________________________
        //this part is if interaction is in oxbash format
        tav::STensor<int,double> ENT(2);  //single particle vector, initilized to zero
        tav::STensor<int,double> INT(6); //two-body
        SM::ReadXMLINT(ENT, INT, xMainNode);
        SM::MakeESPO(ESP,sps,ENT);
        PPInteraction2sp(VSPX, sps, INT);
    }
    int jispfiles=xHamiltonianNode.nChildNode(xsm::JISP);
    if (jispfiles!=0) { //there is jisp interaction
        cerr<<"Interaction files are in JISP format "<<endl;
        XMLNode xJISPNode=xHamiltonianNode.getChildNode(xsm::JISP);
        string JISPFileString;
        ReadXMLData(JISPFileString,"name",xJISPNode);
        int shl;
        ReadXMLData(shl,xsm::maxshell,xJISPNode);
        double hw;
        ReadXMLData(hw,"hw",xJISPNode);
        cerr<<"hw="<<hw<<" maxshell="<<shl<<" file="<<JISPFileString<<endl;
        double massratio=sqrt(938.92/938.093);
        // vector<double> scalings={2*hw/double(N), 0.0, massratio,1.0,  1.0,  1.0};  //requires c++11
        std::vector<double> scalings (6);
        scalings[0]=2*hw/double(N);scalings[1]=0.0; scalings[2]=massratio; scalings[3]=scalings[4]=scalings[5]=1.0;
        //vector<int> ispins={2,                           0 ,  -1,              0, -1,    1   };
        std::vector<int> ispins(6);
        ispins[0]=2; ispins[1]=2; ispins[2]=-1; ispins[3]=0; ispins[4]=-1; ispins[5]=1;
        double tempReadTagDouble = 0;
        
        if (ReadXMLData(tempReadTagDouble,"Trel",xJISPNode)!=1) scalings[0] *= tempReadTagDouble;
        if (ReadXMLData(tempReadTagDouble,"Vrel",xJISPNode)!=1) scalings[1] = tempReadTagDouble;//relative is 0
        if (ReadXMLData(tempReadTagDouble,"Vcoul",xJISPNode)!=1) scalings[2] *= tempReadTagDouble;
        if (ReadXMLData(tempReadTagDouble,"Vpn",xJISPNode)!=1) scalings[3] *= tempReadTagDouble;
        if (ReadXMLData(tempReadTagDouble,"Vpp",xJISPNode)!=1) scalings[4] *= tempReadTagDouble;
        if (ReadXMLData(tempReadTagDouble,"Vnn",xJISPNode)!=1) scalings[5] *= tempReadTagDouble;

        ReadSPINT(VSPX,sps,scalings,ispins,(av::FileName(JISPFileString,"JISPPATH",".INT")).c_str());
    }
    

    if (tagfiles != 0)
    {
        std::cout << "Reading Interactions from TagFiles.\n";
        for (int tagfile_count = 0; tagfile_count < tagfiles; ++tagfile_count)
        {
            XMLNode xTagNode = xHamiltonianNode.getChildNode("tagfile",tagfile_count);
            std::string nameString;
            ReadXMLData(nameString,"fileName",xTagNode);
            
            
            std::cout << "Reading From: " << nameString << std::endl;
//            DataStream.open((av::FileName(sysfile,"SMINTPATH",".INT")).c_str());

            std::cout << "[>] " << av::FileName(nameString.c_str(),"SMINTPATH",".TAG") << std::endl;
            std::ifstream instream(av::FileName(nameString.c_str(),"SMINTPATH",".TAG").c_str());
            
            av::SkipComments(instream);
            long startOfStream = instream.tellg();

            int numberOfEnergies;
            std::string tempString;
            std::vector<double> neut,prot,tot;

            if (av::FileStringPositionFind("Eneutrons:",instream,true) != -1)
            {
                std::cerr << "    [+] Shifting Neutron SPE. " << std::endl;
                instream >> tempString;//save "Eneutrons:"
                instream >> numberOfEnergies;
                neut.resize(numberOfEnergies,0.);
                for (int i = 0; i < numberOfEnergies; ++i)
                    instream >> neut[i] ;
                instream.seekg(startOfStream);//Return the stream where it was to look for other tags.
                tav::STensor<spsint,double>::iterator it = ESP.begin();
                for (; it != ESP.end(); ++it)
                {
                    if( sps[it->first[0]].Tz > 0)
                        it->second += neut[sps[it->first[0]].level];
                    
                }
            }
            
            if (av::FileStringPositionFind("Eprotons:",instream,true) != -1)
            {
                std::cerr << "    [+] Shifting Proton SPE. " << std::endl;
                instream >> tempString;//save "Eprotons:"
                instream >> numberOfEnergies;
                prot.resize(numberOfEnergies,0.);
                for (int i = 0; i < numberOfEnergies; ++i)
                    instream >> prot[i] ;
                instream.seekg(startOfStream);//Return the stream where it was to look for other tags.

                tav::STensor<spsint,double>::iterator it = ESP.begin();
                for (; it != ESP.end(); ++it)
                {
                    if( sps[it->first[0]].Tz < 0)
                        it->second += prot[sps[it->first[0]].level];
                    
                }
            }
            
            if (av::FileStringPositionFind("SPE:",instream,true) != -1)
            {
                std::cerr << "    [+] Shifting SPE. " << std::endl;
                instream >> tempString;//save "SPE:"
                instream >> numberOfEnergies;
                tot.resize(numberOfEnergies,0.);
                for (int i = 0; i < numberOfEnergies; ++i)
                    instream >> tot[i] ;
                instream.seekg(startOfStream);//Return the stream where it was to look for other tags.
                
                tav::STensor<spsint,double>::iterator it = ESP.begin();
                for (; it != ESP.end(); ++it)
                        it->second += tot[sps[it->first[0]].level];
            }
            
//            std::cout << "Extra Neutron SPEs: ";
//            for (int i = 0; i < neut.size(); ++i)
//                std::cout << neut[i] << " ";
//            std::cout << std::endl << "Extra Proton SPEs: ";
//            for (int i = 0; i < prot.size(); ++i)
//                std::cout << prot[i] << " ";
//            std::cout << std::endl;
            
            if (av::FileStringPositionFind("Ecore:",instream,true) != -1)
            {
                std::cerr << "    [+] Adding Core Energy. " << std::endl;
                double CoreEnergy;
                instream >> tempString;
                instream >> CoreEnergy;
                instream.seekg(startOfStream);//Return the stream where it was to look for other tags.

                tav::STensor<spsint,double>::iterator it = ESP.begin();
                for (; it != ESP.end(); ++it)
                {
                    if( sps[it->first[0]].Tz < 0)
                        it->second += CoreEnergy/N;
                    else
                        it->second += CoreEnergy/N;
                }

            }
            if (av::FileStringPositionFind("Particle-Hole Channel",instream,true) != -1)
            {
                std::cerr << "    [+] File is in Particle-Hole Format." << std::endl;
                int linecount = 0,TzType = 0;
                bool selectPN = false;
                std::string tempSTR;
                instream >> tempSTR >> tempSTR;//Read "Particle-Hole Channel"
                if (av::FileStringPositionFind("Type:",instream,true) != -1)
                {
                    selectPN = true;
                    instream >> tempSTR;//Read "Type:"
                    instream >> tempSTR;
                    
                    if (tempSTR == "p" )      {TzType = -1; tempSTR = "proton";}
                    else if (tempSTR == "n" ) {TzType =  1; tempSTR = "neutron";}
                    else if (tempSTR == "pn" ){TzType =  0; tempSTR = "proton-neutron";}
                    else {std::cerr << "    [!] Unknown Type in Tag file... Exiting..." << std::endl;exit(1);}
                    std::cerr << "    [+] Interaction is only in " << tempSTR << " channel." << std::endl;
                }
                if(av::FileStringPositionFind("Lines:",instream,true) != -1)
                {
                    instream >> tempSTR; //Read "Lines: "
                    instream >> linecount;
                    std::cerr << "    [+] " << linecount << " line(s) in file.\n";
                    /////LINES COPIED FROM xsphericalint.cxx/////////
                    /////Only name of XMLNode was changed///////////
                    double myscaling[10];
                    for (int ii=0;ii<10;ii++) myscaling[ii]=1.0;
                    int nscale=xTagNode.nChildNode("scale");
                    for (int s=0; s<nscale; s++) {
                        int rank;
                        getXMLData(rank,"rank", xTagNode.getChildNode("scale",s));
                        double scalerank=atof(xTagNode.getChildNode("scale",s).getText());
                        cerr<<"    [+] Scaling rank "<<rank<<" components with multiplier "<<scalerank<<endl;
                        myscaling[rank]*= scalerank;
                    }
                    /////END COPY-PASTED LINES//////////
                    tav::STensor<int,double> VSPH(6);
                    int* key = new int[6];
                    double mul = 0;
                    for (int ll = 0; ll < linecount; ++ll)
                    {
                        //Read 1234JT
                        for (int xx = 0; xx < VSPH.rank; ++xx)
                        {
                            instream >> key[xx] ;
                            if (xx < 4)//JT don't need reduction.
                                key[xx]--;
                        }
                        //Read value
                        instream >> mul;
                        //Multiply with rank scale and store.
                        VSPH[const_cast<const int*>(key)] = myscaling[2] * mul;
                    }
    //                std::cout << VSPH << "\nxxxxxxxxxxxxxx\n";
                    for (auto it = VSPH.begin(); it != VSPH.end(); ++it)
                        SM::AddQQVspJT(ESP,VSPX,sps,it->first,it->second,selectPN,TzType);
                    
    //                std::cout << ESP << "\n------------\n" << VSPX << std::endl;
                }
            }
            
            

        }
    }
    
    //add J^2 term
    double jadd=0.;
    //funciton returns 1 on fail, so if we do not fail than read
    if (ReadXMLData(jadd,"JJ",xMatrixNode)!=1) {
        cerr<<"Applying JJ term "<<jadd<<endl;
        //my hamiltonian is in ESP and VSPX
        tav::STensor<spsint,double> SPJJ(4);
        tav::STensor<spsint,double> JSP(2);
        //JSP=0.0;
        //SPJJ.clear();
        SM::SPJJ(JSP,SPJJ,sps);
        //SM::SPJJ(ESP,VSPX,sps);
        ESP+=(JSP*=jadd);
        VSPX+=(SPJJ*=jadd);
    }
    
    double tadd=0.;
    //funciton returns 1 on fail, so if we do not fail than read
    if (ReadXMLData(tadd,"TT",xMatrixNode)!=1) {
        cerr<<"Applying TT term "<<tadd<<endl;
        //my hamiltonian is in ESP and VSPX
        tav::STensor<spsint,double> SPJJ(4);
        tav::STensor<spsint,double> JSP(2);
        //JSP=0.0;
        //SPJJ.clear();
        SM::SPTT(JSP,SPJJ,sps);
        //SM::SPJJ(ESP,VSPX,sps);
        ESP+=(JSP*=tadd);
        VSPX+=(SPJJ*=tadd);
    }
    
    //only add lawson term if relevant file is included
    //#ifdef __CM__LAWSON__CXX__
    double beta=0.;
    if (ReadXMLData(beta,"Lawson",xMatrixNode)!=1 || ReadXMLData(beta,"Lawson",xHamiltonianNode)!=1 ) {
        cerr<<"Applying Lawson term with β= "<<beta<<endl;
        
        //Read The Valence Space.
        XMLNode xValenceNode=xMainNode.getChildNode(xsm::valencespace);
        SM::ValenceSpace x;
        ReadXMLValenceSpace(x,xValenceNode);
        
        tav::STensor<int,double> CM(6);
        tav::STensor<spsint,double> CM2(2);
        tav::STensor<spsint,double> VSPCM(4);
        GenerateCMINT(CM,x);
        CM*=beta;
//        MakeCM_LawsonInteraction(sps,CM,beta);
        spsint kk[2];
        for (kk[0]=0;kk[0]<sps.size();kk[0]++)
        {
            kk[1]=kk[0];
            CM2[kk]=beta*(2*sps[int(kk[0])].N+sps[int(kk[0])].L);
        }
        PPInteraction2sp(VSPCM,sps,CM);
        ESP+=(CM2);
        VSPX+=(VSPCM);
    }
    //#endif
    beta=0;
    if (ReadXMLData(beta,"Coulomb",xMatrixNode)!=1 || ReadXMLData(beta,"Coulomb",xHamiltonianNode)!=1)
    {
        
        cerr<<"Applying Coulomb Force with  e^2= "<<beta<<endl;
        //Read The Valence Space.
        XMLNode xValenceNode=xMainNode.getChildNode(xsm::valencespace);
        SM::ValenceSpace x;
        ReadXMLValenceSpace(x,xValenceNode);
        
        tav::STensor<int,double> COULOMB(6);
        tav::STensor<spsint,double> VSPCOULOMB(4);

        GeneraterToTheKINT(COULOMB,x,-1);//Generate Coulomb in spherical Matrix Elements
        COULOMB*=beta;
        PPInteraction2sp(VSPCOULOMB,sps,COULOMB);//Convert to pp format
        //Now set to 0 all p-n elements and all n-n elements. (Proton Tz=-1/2)
        for (tav::STensor<spsint,double>::iterator it=VSPCOULOMB.begin();it!=VSPCOULOMB.end();it++)
            if ((double(sps[it->first[0]].Tz)+double(sps[it->first[1]].Tz)+1)>0.5)
                it->second=0;
//        cout<<"VCOUL:"<<endl<<VSPCOULOMB<<endl;
        VSPX+=(VSPCOULOMB);
    }
    
    beta=0;
    if (ReadXMLData(beta,"Delta",xMatrixNode)!=1 || ReadXMLData(beta,"Delta",xHamiltonianNode)!=1)
    {
        cerr<<"Applying Surface Delta Force with  a= "<<beta<<endl;
        //Read The Valence Space.
        XMLNode xValenceNode=xMainNode.getChildNode(xsm::valencespace);
        SM::ValenceSpace x;
        ReadXMLValenceSpace(x,xValenceNode);
        
        tav::STensor<int,double> DELTA(6);
        tav::STensor<spsint,double> VSPDELTA(4);
        
        GenerateSDIINT(DELTA,x);//Generate Coulomb in spherical Matrix Elements
        DELTA*=beta;
        PPInteraction2sp(VSPDELTA,sps,DELTA);//Convert to pp format
        VSPX+=(VSPDELTA);
    }
    return VSPX.size();
}
#endif //__XREADSPINT__CXX__
