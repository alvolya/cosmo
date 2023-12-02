/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*!
 \author Alexander Volya  <http://www.volya.net>
 \file cosmoxml.cxx
 \ingroup gp_nr
 */
#ifndef __COSMOXML__CXX__
#define __COSMOXML__CXX__
#include <xml/xmlParser.cxx>
using xml::XMLNode;
#include <xsm/schema.h>
#include <sm/spinclass.cxx>
#include <sm/parityclass.cxx>
#include <phy/nuclearnames.cxx>
using std::cout;
using std::cerr;
using std::endl;
#include <xsm/xlabel.cxx>
using SM::LabelSystem;
#include <algorithm>
using  std::replace;
//#include <string>
//using std::replace;
/*Deduce possible parity of the system, if definite, return 0, elese return 1
 */
namespace xsm {
    int SetParity(XMLNode xSystem, int N, XMLNode xValenceSpaceNode) {
        if (xValenceSpaceNode.isEmpty()) FatalError("SetParity, empty ValenceSpace node on input");
        
        int orbits=xValenceSpaceNode.nChildNode(xsm::valencespace_orbital);
        //get first orbit
        int l=0;
        XMLNode xOrbit=xValenceSpaceNode.getChildNode(xsm::valencespace_orbital,0);
        //cout<<xValenceSpaceNode.createXMLString(true);
        try {xml::getXMLData(l,"l", xOrbit);}
        catch (...) {
            //AskParity(xSystem);
            return 1;
        }
        
        //continue with the rest of orbits
        for (int ob=1;ob<orbits;ob++) {
            int ll;
            xOrbit=xValenceSpaceNode.getChildNode(xsm::valencespace_orbital,ob);
            try {xml::getXMLData(ll,"l", xOrbit);}
            catch (...) {
                //AskParity(xSystem);
                return 1;
            }
            //if (SM::parity(l)!=SM::parity(ll)) {
            if ((l&1)!=(ll&1)) {
                //AskParity(xSystem);
                return 1;
            }
        }
        //we got here because all orbits have the same parity l
        if ((l&1)&(N&1)) addAttribute(xSystem,xsm::P, "-"); else addAttribute(xSystem,xsm::P, "+");
        return 0;
    }
    int PHRejections(XMLNode &OutMainNode)
    {
        XMLNode xValenceSpaceNode=OutMainNode.getChildNode(xsm::valencespace);
        if (xValenceSpaceNode.isEmpty()) FatalError("Add rejection, valencespace node is empty");
        XMLNode xSystem=OutMainNode.getChildNode(xsm::system);
        if (xSystem.isEmpty()) xSystem=OutMainNode.addChild(xsm::system);
        
        if(Inquire("Create Particle-Hole Rejection (use with XSHLMBSQRR)",false))
        {
            XMLNode xPHRejection=xSystem.addChild("phrejection");
            int orbits=xValenceSpaceNode.nChildNode(xsm::valencespace_orbital);
            cerr<<"There are 3 types of levels: [0] for empty level, [1] for filled level, [p] for particle level, [h] for hole level."<<endl;
            for (int ob=0;ob<orbits;ob++)
            {
                XMLNode xOrbit=xValenceSpaceNode.getChildNode(xsm::valencespace_orbital,ob);
                
                XMLNode xOrbitNew=xOrbit.deepCopy();
                xPHRejection.addChild(xOrbitNew);
                
                cerr<<"Orbital: "<<xOrbit.getAttribute("name")<<" type= "<<xOrbit.getAttribute("type");
                std::string cinput;
                std::string tyype(xOrbit.getAttribute("type"));
                if (tyype=="pn")
                {
                    cerr<<endl;
                    Inquire(cinput,"Proton (0, 1, p or h): ","0/1/p/h");
//                    updateAttribute(xOrbitNew,"type","p");
                    
                    xOrbitNew.addAttribute("proton",cinput.c_str());
                    cinput.clear();
                    Inquire(cinput,"Neutron (0, 1, p or h): ","0/1/p/h");
                    xOrbitNew.addAttribute("neutron",cinput.c_str());
            
                    

                }
                else if (tyype=="p")
                {
                    Inquire(cinput,"Proton Orbital (0, 1, p or h): ","0/1/p/h");
                    xOrbitNew.addAttribute("proton",cinput.c_str());
                }
                else if (tyype=="n")
                {
                    Inquire(cinput,"Neutron Orbital (0, 1, p or h): ","0/1/p/h");
                    xOrbitNew.addAttribute("neutron",cinput.c_str());
                }
                else
                {
                    cerr<<"Unknown type??? "<<tyype<<std::endl;
                }
            }
            int nmin=0,nmax=0;
            
            Inquire(nmin,"Minimum particle - hole on selected orbits: ");
            Inquire(nmax,"Maximum particle - hole on selected orbits: ");
            addAttribute(xPHRejection,"MINPH", nmin);
            addAttribute(xPHRejection,"MAXPH", nmax);

        }
        return 0;
    }
    int AddRejections(XMLNode xMainNode) {
        XMLNode xValenceSpaceNode=xMainNode.getChildNode(xsm::valencespace);
        if (xValenceSpaceNode.isEmpty()) FatalError("Add rejection, valencespace node is empty");
        XMLNode xSystem=xMainNode.getChildNode(xsm::system);
        if (xSystem.isEmpty()) xSystem=xMainNode.addChild(xsm::system);
        ///here we work on rejectionsc
        while (Inquire("Create rejection (to be processed by Xsysmbs)",false)) {
            XMLNode xRejection=xSystem.addChild("rejection");
            int orbits=xValenceSpaceNode.nChildNode(xsm::valencespace_orbital);
            XMLNode xOrbit;
            int NMAX, NMIN;
            NMAX=NMIN=0;
            for (int ob=0;ob<orbits;ob++) {
                xOrbit=xValenceSpaceNode.getChildNode(xsm::valencespace_orbital,ob);
                cerr<<"Orbital: "<<xOrbit.getAttribute("name")<<" type= "<<xOrbit.getAttribute("type");
                std::string cinput;
                if (Inquire("include ",false)) {
                    //XMLNode xOrbitNew=XMLNode::parseString(xOrbit.createXMLString(false));
                    XMLNode xOrbitNew=xOrbit.deepCopy();
                    xRejection.addChild(xOrbitNew);
                    Inquire(cinput,"neutron, proton, or both (n, p, or pn): ","pn");
                    //AskUser("neutron, proton, or both (n, p, or pn): ", cinput);
                    //xOrbitNew.deleteAttribute("type");
                    //xOrbitNew.addAttribute("type",cinput.c_str());
                    updateAttribute(xOrbitNew,"type",cinput);
                    //-for fun count maximum number of particles
                    SM::spin ospin;
                    int nn=int(getXMLData(ospin,xsm::J,xOrbitNew))+1;
                    if (cinput=="pn") nn*=2;
                    NMAX+=nn;
                } //end of add
            } //end looping over orbits
            //AskUser("Minimum N on selected orbits: ", NMIN);
            std::cout<<"The number of nucleons in selected orbits is in ["<<NMIN<<","<<NMAX<<"]"<<std::endl;
            Inquire(NMIN,"Minimum N on selected orbits: ");
            addAttribute(xRejection,"NMIN", NMIN);
            Inquire(NMAX,"Maximum N on selected orbits: ");
            //AskUser("Maximum N on selected orbits: ", NMAX);
            addAttribute(xRejection,"NMAX", NMAX);
            //AskUser("Add another rejection (y/n) ", go);
        }
        return 0;
    } //end of function
    
    int Rejections(XMLNode xMainNode) {
        XMLNode xHamiltonian=xMainNode.getChildNode(xsm::hamiltonian);
        if (xHamiltonian.isEmpty()) CriticalError("Rejection, hamiltonian is empty");
        int nrej=xHamiltonian.nChildNode("rejection");
        if (nrej==0) {return AddRejections(xMainNode);}
        //std::string cinput;
        //AskUser("do you want to use rejections provided with the Hamiltonian (y/n): ", cinput);
        if (Inquire("do you want to use rejections provided with the Hamiltonian: ",false))    {
            //copy all rejections from the Hamiltonian
            XMLNode xSystem=xMainNode.getChildNode(xsm::system);
            if (xSystem.isEmpty()) xSystem=xMainNode.addChild(xsm::system);
            for (int r=0;r<nrej;r++) {
                //xSystem.addChild(xHamiltonian.getChildNode("rejection")); //this will move (always first)
                xSystem.addChild(xHamiltonian.getChildNode("rejection",r).deepCopy()); //this will copy
                //now we need to copy
            }
            //add more rejections if needed
        } //end if to see if one wants hamiltonian rejections
        return AddRejections(xMainNode);
    }

int HWRejections(XMLNode xSystem) {
if(Inquire("Create hw rejection (use with Xsysmbs)",false))
           {
               XMLNode xHWRejection=xSystem.addChild("hwrejection");
               int HWMAX=0;
               int HWMIN=0;
               Inquire(HWMAX,"Maximum number of oscillator qunta above minimum: ");
               Inquire(HWMIN,"Minimum number of oscillator qunta above minimum: ");
               addAttribute(xHWRejection,"HWMAX", HWMAX);
               addAttribute(xHWRejection,"HWMIN", HWMIN);
           }
    return 0;
}
    
    
    int SelectModelSpace(
                         XMLNode &xNodeSP, ///<output node containing model space
                         XMLNode xMainNode, ///<input main int.xml node
                         int N, ///<neutron number
                         int Z ///<proton number
    )
    {
        char *defaultmodelsp;
        int NMIN, NMAX;
        int ZMIN, ZMAX;
        int nsp=xMainNode.nChildNode("modelspace");
        
        if (nsp!=0) { //this loop is only done for a big database search, if there are no model spaces then xNodeSP=xMainNode
            cerr<<"number of spaces "<<nsp<<endl;
            for (int i=0;i<nsp;i++) {
                //NMIN=xMainNode.getChildNode("modelspace/valencespace").getAttributeValue()
                xNodeSP = xMainNode.getChildNode("modelspace",i);
                XMLNode xNode2 = xNodeSP.getChildNodeByPath("valencespace/core");
                //cerr<<xNode2.getAttribute(xsm::nucleus_N)<<endl;
                getXMLData(NMIN,xsm::nucleus_N,xNode2);
                getXMLData(ZMIN,xsm::nucleus_Z,xNode2);
                //NMIN=atoi(xNode2.getAttribute(xsm::nucleus_N));
                //ZMIN=atoi(xNode2.getAttribute(xsm::nucleus_Z));
                string AMIN = xNode2.getAttribute("name");
                xNode2 = xNodeSP.getChildNodeByPath("valencespace/outercore");
                getXMLData(NMAX,xsm::nucleus_N,xNode2);
                getXMLData(ZMAX,xsm::nucleus_Z,xNode2);
                //cerr<<xNode2.getAttribute(xsm::nucleus_N)<<endl;
                //NMAX=atoi(xNode2.getAttribute(xsm::nucleus_N));
                //ZMAX=atoi(xNode2.getAttribute(xsm::nucleus_Z));
                string AMAX = xNode2.getAttribute("name");
                //cerr<<NMAX<<" "<<ZMAX<<" "<<NMIN<<" "<<ZMIN<<endl;
                if ((N<=NMAX)&&(Z<=ZMAX)&&(N>=NMIN)&&(Z>=ZMIN)) {
                    cout<<"you can use modelspace "<<xNodeSP.getChildNode(xsm::valencespace).getAttribute("name")<<"  \t "<<AMIN<<" - "<<AMAX<<endl;
                } //end if
            } //end sp loop
            
            
            /*
             string modelsp;
             AskUser("select model space ", modelsp);
             //with this loop we set the desired model space
             for (int i=0;i<nsp;i++) {
             xNodeSP=xMainNode.getChildNode("modelspace",i);
             string tmpmodelsp=xNodeSP.getChildNode(xsm::valencespace).getAttribute("name");
             //cout<<tmpmodelsp<<endl;
             if (tmpmodelsp==modelsp) break;
             }
             
             */
            char modelsp[256];
            bool correctSpace = false;
            while (!correctSpace)
            {
                AskUser("select model space ", modelsp);
                //search for hamiltonian with this name
                xNodeSP=xMainNode.getChildNodeWithAttribute("modelspace","name",modelsp);
//                if (xNodeSP.isEmpty()) CriticalError("There is no such model space");
                if (xNodeSP.isEmpty()) std::cerr << "There is no such model space." << std::endl;
                else correctSpace = true;
            }
        } //end search for model spaces
        else {xNodeSP=xMainNode;}
        return 0;
    }
    
    int SelectModelHamiltonian(
                               XMLNode &xNodeINT, ///<output main hamiltonian node
                               XMLNode xNodeSP ///<input node containing model space
    ) {
        int nint=xNodeSP.nChildNode(xsm::hamiltonian);
        cerr<<"Available Hamiltonians"<<endl;
        for (int i=0;i<nint;i++) {
            //cerr<<i<<endl;
            cout<<"--------"<<xNodeSP.getChildNode(xsm::hamiltonian,i).getAttribute("name")<<"--------"<<endl;
            if (!xNodeSP.getChildNode(xsm::hamiltonian,i).getChildNodeByPath("reference").isEmpty())
                cout<<xNodeSP.getChildNode(xsm::hamiltonian,i).getChildNodeByPath("reference").getText()<<endl;
            if (!xNodeSP.getChildNode(xsm::hamiltonian,i).getChildNode(xsm::hamiltonian_comment).isEmpty())
                cout<<xNodeSP.getChildNode(xsm::hamiltonian,i).getChildNode(xsm::hamiltonian_comment).getText()<<endl<<endl;
        }
        
        char modelham[256];
        
        if (nint==1) xNodeINT=xNodeSP.getChildNode(xsm::hamiltonian);
        else {
            bool correctHamiltonian = false;
            while (!correctHamiltonian)
            {
                AskUser("select Hamiltonian ", modelham);
                //search for hamiltonian with this name
                xNodeINT=xNodeSP.getChildNodeWithAttribute(xsm::hamiltonian,"name",modelham);
//                if (xNodeINT.isEmpty()) CriticalError("There is no such hamiltonian");
                if (xNodeINT.isEmpty()) std::cerr << "There is no such hamiltonian" << std::endl;
                else correctHamiltonian = true;
            }
        }
      return 0;
    }
    
    
    //Make an xml header with a system and with Hamiltonian
    int CoSMoXMLHeader(
                       XMLNode &xOutMainNodeX //return node with header
    ) {
        xOutMainNodeX=XMLNode::createXMLTopNode("xml",TRUE);
        addAttribute(xOutMainNodeX,"version","1.0");
        XMLNode xOutMainNode=xOutMainNodeX.addChild("cosmo");
        addAttribute(xOutMainNode,"version","4.1");
#ifdef COSMO_REVISION
        addAttribute(xOutMainNode,"revision",COSMO_REVISION);
#endif
        string compiledate=__DATE__;
        replace(compiledate.begin(), compiledate.end(), ' ', '/');
        //replace(compiledate.begin(), compiledate.end(), ' ', '/');
        addAttribute(xOutMainNode,"date", compiledate.c_str());
        
        return 0;
    }
    
    //return name as a string
    string AddXMLSystem(
                        XMLNode &xOutMainNode, //cosmo main
                        int N,
                        int Z
                        )
    {
        
        XMLNode  xNode2=xOutMainNode.getChildNodeByPath("valencespace/core");
        //cerr<<xNode2.getAttribute(xsm::nucleus_N)<<endl;
        int NMIN;
        int ZMIN;
        getXMLData(NMIN,xsm::nucleus_N,xNode2);
        getXMLData(ZMIN,xsm::nucleus_Z,xNode2);
        int AV=N+Z-NMIN-ZMIN;
        XMLNode xSystem=xOutMainNode.addChild(xsm::system);
        //XMLNode x2System=xSystem.addChild("parameters");
        addAttribute(xSystem,xsm::valence_A,AV);
        addAttribute(xSystem,xsm::valence_N,N-NMIN);
        addAttribute(xSystem,xsm::valence_Z,Z-ZMIN);
        updateAttribute(xSystem,xsm::nucleus_A,N+Z);
        updateAttribute(xSystem,xsm::Tz,SM::spin(N-Z));
        SM::spin Jz=SM::spin(AV%2);
        InquireAttribute(Jz,xsm::Jz, "Enter spin projection Jz", xSystem);
        if(Inquire("Would you like to shift hamiltonian matrix by J+J- term",false))
        {
            XMLNode xMatrixNode=xSystem.addChild(xsm::system_matrix);
            double jadd=0.0;
            tav::Inquire(jadd,"Enter scaling for J+J- part ");
            addAttribute(xMatrixNode,"JJ",jadd);
        }
        //if (SetParity(xSystem,N+Z-NMIN-ZMIN,xOutMainNode.getChildNode(xsm::valencespace))) AskParity(xSystem);
        //this sets parity if possible
        //cout<<xOutMainNode.createXMLString(true);
        if (SetParity(xSystem,N+Z-NMIN-ZMIN,xOutMainNode.getChildNode(xsm::valencespace))){
            string parity;
            AskAttribute(parity,xsm::P, "Enter parity (+ or -): ", xSystem);
        };
        //char parity[2];
        //AskUser("Enter parity (+ or -) ", parity);
        //x2System.addAttribute(xsm::P, parity);
        HWRejections(xSystem);
        PHRejections(xOutMainNode);
        Rejections(xOutMainNode);
        //overwrite argument
        string tmpname;
        XMLNode xNodeINT=xOutMainNode.getChildNode(xsm::hamiltonian);
        ReadXMLData(tmpname,"name",xNodeINT);
        stringstream nuclearname;
        nuclearname<<N+Z<<phy::NuclearSymbol[Z];
        string nucintname=nuclearname.str()+"_"+tmpname;
        string tmpname2=tmpname=nucintname+"_"+LabelSystem(xSystem);
        Inquire(tmpname, "Name the sytem as ");
        addAttribute(xSystem,"name",tmpname.c_str());
        if (tmpname2!=tmpname) return tmpname; //test if there is a change in name
        return nucintname;
    }

    string AddXMLSystem(
                        XMLNode &xOutMainNode, //cosmo main
                        int N,
                        int Z,
                        SM::spin Jz,
                        double jadd
                        )
    {
        
        XMLNode  xNode2=xOutMainNode.getChildNodeByPath("valencespace/core");
        //cerr<<xNode2.getAttribute(xsm::nucleus_N)<<endl;
        int NMIN;
        int ZMIN;
        getXMLData(NMIN,xsm::nucleus_N,xNode2);
        getXMLData(ZMIN,xsm::nucleus_Z,xNode2);
        int AV=N+Z-NMIN-ZMIN;
        XMLNode xSystem=xOutMainNode.addChild(xsm::system);
        //XMLNode x2System=xSystem.addChild("parameters");
        addAttribute(xSystem,xsm::valence_A,AV);
        addAttribute(xSystem,xsm::valence_N,N-NMIN);
        addAttribute(xSystem,xsm::valence_Z,Z-ZMIN);
        updateAttribute(xSystem,xsm::nucleus_A,N+Z);
        updateAttribute(xSystem,xsm::Tz,SM::spin(N-Z)); 
        addAttribute(xSystem,xsm::Jz, Jz);
        XMLNode xMatrixNode=xSystem.addChild(xsm::system_matrix);
        addAttribute(xMatrixNode,"JJ",jadd);
        
        //if (SetParity(xSystem,N+Z-NMIN-ZMIN,xOutMainNode.getChildNode(xsm::valencespace))) AskParity(xSystem);
        //this sets parity if possible
        //cout<<xOutMainNode.createXMLString(true);
        if (SetParity(xSystem,N+Z-NMIN-ZMIN,xOutMainNode.getChildNode(xsm::valencespace))){
            string parity;
            AskAttribute(parity,xsm::P, "Enter parity (+ or -): ", xSystem);
        };
        //char parity[2];
        //AskUser("Enter parity (+ or -) ", parity);
        //x2System.addAttribute(xsm::P, parity);
        Rejections(xOutMainNode);
        //overwrite argument
        string tmpname;
        XMLNode xNodeINT=xOutMainNode.getChildNode(xsm::hamiltonian);
        ReadXMLData(tmpname,"name",xNodeINT);
        stringstream nuclearname;
        nuclearname<<N+Z<<phy::NuclearSymbol[Z];
        string nucintname=nuclearname.str()+"_"+tmpname;
        tmpname=nucintname+"_"+LabelSystem(xSystem); 
        Inquire(tmpname, "Name the sytem as ");
        addAttribute(xSystem,"name",tmpname.c_str());
        return nucintname;
    }
    
} //end namespace xsm

#endif //__COSMOXML__CXX__

