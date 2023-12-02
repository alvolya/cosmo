/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/**
 
 \author Alexander Volya  <http://www.volya.net>
 \file Xsysmbs.cpp
 \brief Create many-body states by partitioning, allow rejections
 \ingroup gp_mbs
 \todo
 - Separate partitions into separate class
 - Make2 with partitions has a problem
 - Assembling partitions for all levels into a state is done in a primitive way by using vector access and iterating over all possibilities.
 */
/*
 08/07/08 rejections are included in merging partitions, prevent memory problem
 09/27/08 status and wall times
 */
//#define STATUS yes
//#define WALLTIME 300
//#include <sm/ValenceSpace.cxx>
#if __cplusplus <= 199711L
typedef unsigned short spsint; //this part is for old c++98
#else
#include <cstdint>
typedef uint16_t spsint;
#endif
typedef unsigned long long  uint_nbasis;
#include <xsm/xSingleParticleStates.cxx>
#include <iostream>
typedef unsigned char uschar;
#include "PartitionClass.cxx"
#include "SpinPartitionClass.cxx"
//#include <tav/sort.cxx>
#include <tav/stensor.cxx>
#include <sm/label.cxx>
#include <inquire.cxx>
using tav::Inquire;

//#include <status.h>
using SM::AddChar;
using namespace std;

#include <av/argparser.cxx>



/*-----------MAIN---------------------------------------*/
int main (int argc, char** argv){
av::argparser args(argc,argv);
args.SetDescription("Create many-body states by partitioning, allow rejections.");
args.SetDefault("XML file containing a parent system description.",2,2);
args.SetFlag("o","output","Output XML file, if provided the default input is unchanged.",0,2);
if(args.ParseArguments()) return 1;
  //unsigned char dfd;
  /*----------Check that interaction file is fine----------*/
  string sysfile=args(1);
 
 
 //make system
XMLNode xMainNode=XMLNode::openFileHelper(sysfile.c_str(), xsm::cosmo); //we start by opening xml input file
 SM::SingleParticleStates sps; Make(sps,xMainNode.getChildNode(xsm::valencespace));
    
    // here we use single particle states, it is more than needed but certain to provide compatibility.
    std::map<int,int> spinsprotons;
    std::map<int,int> spinsneutrons;
    std::map<int,int> pxmlindexmap; //this is used to map an xml index to a level protons
    std::map<int,int> nxmlindexmap; //this is used to map an xml index to a level neutrons
   // std::map<int,int>::iterator i_sp;
    //store last sps with these qn for protons and neutrons so spinsneutrons[level]-> last id in sps
    for (int is=0;is<sps.size();is++) {
        if (sps[is].Tz>=0) {spinsneutrons[sps[is].level]=is;
            //cerr<<"neutron "<<is<<" "<<sps[is].level<<endl;
        }
        else spinsprotons[sps[is].level]=is;
    }
//    cerr<<"we filled spinsneutrons "<<endl;
//    for (auto ii=spinsneutrons.begin(); ii!=spinsneutrons.end(); ++ii) cerr<<ii->first<<" "<<ii->second<<endl;
//    cerr<<"we filled spinsprotons "<<endl;
//    for (auto ii=spinsprotons.begin(); ii!=spinsprotons.end(); ++ii) cerr<<ii->first<<" "<<ii->first<<endl;

    
    //create array for spins
    int levels=spinsneutrons.size()+spinsprotons.size();
    int *spins=new int [levels];
    //* Parity rejection -----------------------------
    bool *nparity= new bool [levels];
    vector<int> quanta(levels);
    vector<int> nquanta; //quanta for neutrons
    vector<int> pquanta; //quanta for protons
    //   for (int l=0;l<x.size();l++) {nparity[l]=(x[l].L&1); nparity[l+x.size()]=(x[l].L&1);}
    // cout<<(x[l].P&1)<<" \n";
    // cerr<<x[l].L<<endl;
    //fill array for spins
    //cerr<<spinsneutrons.size()<<endl;
    int ll=0; //just incremental index filling all levels with spin and parity
    for (auto i_sp=spinsneutrons.begin();i_sp!=spinsneutrons.end();i_sp++) {
        spins[ll]=sps[i_sp->second].J+1;
        nparity[ll]=sps[i_sp->second].P;
        quanta[ll]=2*sps[i_sp->second].N+sps[i_sp->second].L;
        nquanta.push_back(2*sps[i_sp->second].N+sps[i_sp->second].L);
        nxmlindexmap[sps[i_sp->second].level]=ll;
        ll++;
        //cerr<<i_sp->first<<" "<<i_sp->second<<" "<<sps[i_sp->second].level<<endl;
    }
    //cerr<<spinsprotons.size()<<endl;
    for (auto i_sp=spinsprotons.begin();i_sp!=spinsprotons.end();i_sp++) {
        spins[ll]=sps[i_sp->second].J+1;
        nparity[ll]=sps[i_sp->second].P;
        quanta[ll]=2*sps[i_sp->second].N+sps[i_sp->second].L;
        pquanta.push_back(2*sps[i_sp->second].N+sps[i_sp->second].L);
        pxmlindexmap[sps[i_sp->second].level]=ll;
        ll++;
        //cerr<<i_sp->first<<" "<<i_sp->second<<" "<<sps[i_sp->second].level<<endl;
    }
    
    //for (i_sp=nxmlindexmap.begin();i_sp!=nxmlindexmap.end();i_sp++) {
    //cerr<<i_sp->first<<" "<<i_sp->second<<" "<<endl;
    //}
    
    //return 1;
    //for (int l=0;l<levels;l++) {cerr<<quanta[l]<<endl;}
    //return 1;
    //spins[0]=2; spins[1]=4; spins[2]=6; //USD shell model
    //spins[3]=2; spins[4]=4; spins[5]=6;
    /*
     Input to the following part of the code
     spins[levels] array of all 2j+1
     the first levels/2 elements are for neutrons
     the second levels/2 elements are for protons
     */
    
    //for convenience we create an offset array.
    int *offset=new int [levels];
    offset[0]=0;
    for (int l=1;l<levels;l++) offset[l]=offset[l-1]+spins[l-1];
    
    
    
    
    SM::QN x;
    
    //data on number of particles
    XMLNode xSystemNode=xml::SelectNode(xMainNode,xsm::system,"name");
    //xMainNode.getChildNode(xsm::system, -1);
    char systemname[256];
    AskAttribute(systemname,"name","Enter system name: ",xSystemNode);
    cerr<<"Processing system : "<<systemname<<endl;
    AskAttribute(x.N,xsm::valence_A, "Enter valence particle number N:", xSystemNode);
    x.Jz=SM::spin(x.N%2); x.Tz=SM::spin(x.N%2);
    int Z;
    int N;
    InquireAttribute(x.Jz,xsm::Jz, "Enter spin projection Jz", xSystemNode);
    try {
        getXMLData(N,xsm::valence_N,xSystemNode);
        getXMLData(Z,xsm::valence_Z,xSystemNode);
    }
    catch (...) {
        //if we use old xml then this data should be avalable
        InquireAttribute(x.Tz,xsm::Tz, "Enter isospin projection Tz=(N-Z)/2 in valence space ", xSystemNode);
        Z=(x.N-x.Tz)/2;
        N=(x.N+x.Tz)/2;;
    }
    cerr<<"Valence particles N="<<N<<" Z="<<Z<<endl;
    int A=N+Z;
    
    char parity[2];
    AskAttribute(parity,xsm::P, "Enter parity (+ or -): ", xSystemNode);
    x.P=(parity[0]=='-');
    
    
    
    
    //char filename[256];
    uint_nbasis stotal=0;
    
    //make an initial distribution
    //  if (sps.size()<A) {cerr<<"fermions do not fit in this space"<<endl; return 0;}
    
    spsint *mystate =new  spsint [A];
    tav::STensor<spsint, char> basis(A);

    Partition NX(N,spinsneutrons.size(),spins); //make neutron partition

    //  cerr<<NX.size()<<endl;
    // Pause();
  
    Partition PX(Z,spinsprotons.size(),spins+spinsneutrons.size()); //make proton partition
  
    //cerr<<PX.size()<<endl;
    //Pause();
    int M=x.Jz; //2*spin projection
    //Partition PP(N,levels,spins);
    // Partition PP(NX,PX);
    Partition PP;
    //-------------Rejection file -------------------
    //ShowStatus();
    XMLNode xManyBody=xSystemNode.getChildNode(xsm::system_basisstates,-1);
    if (!(xManyBody.isEmpty()))
    {
        if(!(Inquire("There are already manybody states associated with this system. Overwrite?",false)))
        {
            exit(0);
        }
    }
    if (xManyBody.isEmpty()) xManyBody=xSystemNode.addChild(xsm::system_basisstates,-1);
    
    int Ncore;
    int Pcore;
    //determine rejection by the number of qunta
    int nhwrej=xSystemNode.nChildNode("hwrejection");
    if (nhwrej>1) cerr<<"Only one rejection for hw is used "<<endl;
    int hwmax=0;
    int hwmin=0;
    bool bhwmax=false;
    try {
        getXMLData(hwmax,"hwmax",xSystemNode);
        hwmin=0;
        bhwmax=true; //set hwmax
        cerr<<"Process hw rejection "<<hwmax<<" "<<hwmin<<endl;
    }
    catch (...) {
        bhwmax=false;
    }
    if (nhwrej>0) {
    xml::getXMLData(hwmax,"HWMAX", xSystemNode.getChildNode("hwrejection"));
    xml::getXMLData(hwmin,"HWMIN", xSystemNode.getChildNode("hwrejection"));
    cerr<<"Process hw rejection "<<hwmax<<" "<<hwmin<<endl;
    }
    if ((nhwrej>0)||bhwmax) {
        Ncore=CoreQuanta(N,+1,sps);
        Pcore=CoreQuanta(Z,-1,sps);
    cerr<<"Core neutron quanta "<<Ncore<<" Core proton quanta "<<Pcore<<endl;
        if ((Ncore<0)||(Pcore<0)) {cerr<<"A problem with quanta"; return 1;}
        NX.MakeQ(Ncore,Ncore+hwmax,nquanta);
        PX.MakeQ(Pcore,Pcore+hwmax,pquanta);
    }
    else {
        NX.Make();
        PX.Make();
    }//selection of hwrejections are defined
    
    int nrej=xSystemNode.nChildNode("rejection");
    cerr<<"The number of rejections to process is "<<nrej<<endl;
    if (nrej!=0) {
        std::vector< std::vector<bool> >rej;
        std::vector<bool> currentrejection(levels);
        std::vector<int> rmin(nrej);
        std::vector<int> rmax(nrej);
        for (int i=0;i<nrej;i++){
            XMLNode xRejection=xSystemNode.getChildNode("rejection", i);
            //rmin[i]=xRejection.getAttribute("NMIN");
            xml::getXMLData(rmin[i],"NMIN", xRejection);
            xml::getXMLData(rmax[i],"NMAX", xRejection);
            //rmax[i]=xRejection.getAttribute("NMAX");
            for (int l=0;l<levels;l++) currentrejection[l]=false; //initilize
            int rob=xRejection.nChildNode(xsm::valencespace_orbital);
            //loop over orbitals included in rejection
            cerr<<"Rejection: "<<i<<" NMIN="<<rmin[i]<<"  NMAX="<<rmax[i]<<" orbitals involved= "<<rob<<endl;
            for (int rr=0;rr<rob;rr++) {
                XMLNode xOrbit=xRejection.getChildNode(xsm::valencespace_orbital, rr);
                int ix;
                xml::getXMLData(ix,"index", xOrbit);
                // const char *t=xOrbit.getAttribute("type");
                ix--; //old oxbash compatibility //note -1 for c-format
                if (xOrbit.getAttribute("type")!=NULL) {
                    std::string type=xOrbit.getAttribute("type");
                    //cerr<<"here"<<endl;
                    if (type=="pn") {currentrejection[pxmlindexmap[ix]]=true; currentrejection[nxmlindexmap[ix]]=true;}
                    if (type=="p") {currentrejection[pxmlindexmap[ix]]=true;}
                    if (type=="n") {currentrejection[nxmlindexmap[ix]]=true;}
                }
                else
                    currentrejection[nxmlindexmap[ix]]=true;
                
                // for (int l=0;l<levels;l++) rejfile>>rej[i][l];
            }
            rej.push_back(currentrejection);
            //rejfile>>rmin[i];
            //rejfile>>rmax[i];
        } //go over rejections
        
        
        MakePartition(PP, NX, PX, rej, rmin, rmax);
    }
    else  MakePartition(PP, NX, PX);
    
    
    
    //ShowStatus();
    // cout<<PP<<endl;
    // Pause();
    
    
    ParityPurge(PP,nparity,x.P);

//final oscillator quanta  purge
   if ((nhwrej>0)||bhwmax) QuantaMinMaxPurge(PP,quanta, Ncore+Pcore+hwmin,Ncore+Pcore+hwmax);
                 
    
    
    
    
    
     std::cout<<"There are "<<PP.size()<<" partitions"<<std::endl;
    //cout<<PP.Make()<<endl;
    //cerr<<"Here is a list of partitions (particles on each level)"<<endl;
    //cerr<<PP<<endl;
    //cerr<<"---------------------"<<endl;
    //for each partition in PP we need to make m-partitions
    //int stotal=0;
    int *mmax=new int [levels];
    int *mmin=new int [levels]; //let us add minimum spin
    int itCounter=0;
    for (Partition::iterator ip=PP.begin(); ip!=PP.end(); ip++) {
        int MM=M;
        for (int i=0;i<levels;i++)
        {
            //mmax[i]=((*ip)[i]*(2*(spins[i])-(*ip)[i]-1))/2;
            /* Maximum \tilde{M}= N\Omega -N(N+1)/2, this is shifted value
             actual max m is N(2j+1-N)/2, this unshifted by jn
             note spins=2j+1
             we cacluate below as plus minus
             */
            mmax[i]=((*ip)[i]*(spins[i]-(*ip)[i])+(*ip)[i]*(spins[i]-1))/2;
            mmin[i]=(-(*ip)[i]*(spins[i]-(*ip)[i])+(*ip)[i]*(spins[i]-1))/2;
            MM+=(*ip)[i]*(spins[i]-1); //we count spins not from zero but as MM-> M+j
        }
        MM/=2;
        //cerr<<"Partitions M-scheme M="<<MM<<endl;
        //for (int i=0;i<levels;i++) cerr<<spins[i]<<" "<<int((*ip)[i])<<" "<<mmax[i]<<" "<<mmin[i]<<endl;
        //cerr<<endl;
        Partition MP(MM,levels,mmax,mmin);///make M projection partition
        MP.Make();
        // cerr<<MP<<endl;
        // Pause();
        for (Partition::reverse_iterator im=MP.rbegin(); im!=MP.rend(); ++im) {
            /* Given data on particle number and spins, generate all data*/
            int mstates=1;
            SpinPartition *SSP= new SpinPartition [levels];
            for (int l=0;l<levels;l++) {
                SSP[l].Make((*im)[l],(*ip)[l], spins[l]);
            }
            //------------Here We have all states
            for (int l=0;l<levels;l++) mstates*=SSP[l].size();
            for (int ii=0;ii<mstates;ii++) {
                int itmp=ii;
                int nnp=0;
                for (int l=0;l<levels;l++) {
                    for (int pp=0;pp<(*ip)[l];pp++)
                    {  mystate[nnp]=SSP[l][itmp%SSP[l].size()][pp]+offset[l]; nnp++; }
                    itmp/=SSP[l].size();
                    //cout<<offset[l]<<" ";
                }
                //cout<<endl;
                //   tav::QuickSort(A,mystate); //sort partition
                /*
                 outf.write((char*)(mystate), A*sizeof(spsint));
                 */
                /*
                 for (int pp=0;pp<A;pp++) cout<<int(mystate[pp])<<" ";
                 cout<<" : ";
                 for (int pp=0;pp<levels;pp++) cout<<(*ip)[pp]<<" ";
                 cout<<" : ";
                 for (int pp=0;pp<levels;pp++) cout<<(*im)[pp]<<" ";
                 cout<<endl;
                 */
                //ShowStatus();
                basis[const_cast<const spsint*>(mystate)]=1;
            }
            //-------------Finish with all states
            delete [] SSP;
            stotal+=mstates;
            //if (basis.size()>100000000) FatalError("The number of m-scheme states limited to 100M");
        }
        
        std::cout << "\rDone: " << itCounter << "/" << PP.size() << "\t(" << double(itCounter)*100./PP.size() << "%) MBS: " << basis.size() << "            "<< std::flush;
        itCounter++;
    }
    //cout<<stotal<<endl;
    //  cout<<x.tag.c_str()<<endl;
    cout<<"number of mbs "<<basis.size()<<endl;
    //now I need to overwrite n in file
    ofstream outf;
    outf.open(AddChar(systemname,".mbs"), ios::binary);
    outf.write((char*)(&A),sizeof(int));
    outf.write((char*)(&stotal),sizeof(uint_nbasis)); //write zero for now
    //  outf.seekp(0); //go to the beguinning of the file
    //  outf.write((char*)(&A),sizeof(int));
    //  outf.write((char*)(&stotal),sizeof(uint_nbasis)); //overide new n and x.N
    for (tav::STensor<spsint, char>::iterator ii=basis.begin(); ii!=basis.end(); ii++) 
        outf.write((char*)(ii->first), A*sizeof(spsint));
    outf.close();
    if (xManyBody.getAttribute("name")!=NULL) xManyBody.deleteAttribute("name");
    if (xManyBody.getAttribute("size")!=NULL) xManyBody.deleteAttribute("size");
    if (xManyBody.getAttribute("creator")!=NULL) xManyBody.deleteAttribute("creator");

    addAttribute(xManyBody,"name",AddChar(systemname,".mbs"));
    addAttribute(xManyBody,"size",int(basis.size()));
    addAttribute(xManyBody,"creator","Xsysmbs");
    // cout<<basis.size()<<endl;
    delete [] mystate;
    delete [] mmax;
    delete [] mmin;
    delete [] spins;
    if (args.isSet("o",1)) xMainNode.writeToFile(args("o",1).c_str());
else xMainNode.writeToFile(sysfile.c_str());

    return 0;
}

