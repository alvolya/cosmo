/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file XSHLEMB.cpp
\brief Electric, magnetic, Gamow-Teller transition rates, includes sumrules
\ref ElectricSPMultipole()
\todo Check magnetic projection sign convetion of creation/annihilation operator
*/
/*
12/21/12 BGT tested and is OK
*/

//#include <cstdint>
#include <cmath>
#if __cplusplus <= 199711L
typedef unsigned short spsint;
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
#include <xsm/xManyBodyStates.cxx>
#include <xsm/xsphericalint.cxx>
#include <sm/MakeSPO.cxx>
#include <xsm/xState.cxx>
#include <sm/EMspherical2sp.cxx>
#include <tav/ThreeJSymbolM.cxx>
using tav::ThreeJSymbol;
#include <sm/PPMatrix.cxx>
//#include <sm/supplemental.cxx>
#include <xsm/xmlemb.cxx>


//Taken from Ring Schuck Appendix B p 592.
double WeisskopfEstimate(XMLNode &xoper,int A, double CoSMoValue)
{
    //First Read E/M transition from xoper.
    string operatorname(xoper.getAttribute("name"));
    char type=operatorname[0];

    if (type!='M' && type!='E') return -1.0;//Weisskopf only for E/M transitions
    int L=atoi(xoper.getAttribute("J"));//multipolarity
    
    if (type=='E')
    {
        double hw = 45. * std::pow(double(A),-1./3.) - 25. * std::pow(double(A),-2./3.);
        
//        double b=1.01*1.01*pow(double(A),1./3.);//OscillatorLength squared
        double b = phy::HBARC* phy::HBARC / ( hw * 938.9186945); //mass defined as avverage of proton neutron.
        double ws=pow(1.2,2*L)*(3./(L+3.))*(3./(L+3.))*pow(double(A),2*L/3.)/(4*M_PI);//Weisskopf Estimate
        return pow(b,L)*CoSMoValue/ws;//convert to (e fm)^2L and divide by Weisskopf Estimate.
    }
    else if (type=='M')
    {
        double hw = 45. * std::pow(double(A),-1./3.) - 25. * std::pow(double(A),-2./3.);
//        double b=1.01*1.01*pow(double(A),1./3.);//OscillatorLength squared
        double b = phy::HBARC* phy::HBARC /( hw * 938.9186945); //mass defined as avverage of proton neutron.
        double ws=10.*pow(1.2,2*L-2)*(3./(L+3.))*(3./(L+3.))*pow(double(A),(2*L-2)/3.)/M_PI;
        return pow(b,L-1)*CoSMoValue/ws;
    }
    //If this point is reached something went horribly wrong.
    cerr<<"Something is wrong in Weisskopf routine"<<endl;
    return -1;
}
//Taken from simple math.
double BLUnitConversion(XMLNode &xoper,int A, double CoSMoValue)
{
    //First Read E/M transition from xoper.
    string operatorname(xoper.getAttribute("name"));
    char type=operatorname[0];
    
    if (type!='M' && type!='E') return -1.0;//Works only for E/M transitions
    int L=atoi(xoper.getAttribute("J"));//multipolarity
    
    if (type=='E')
    {
        double hw = 45. * std::pow(double(A),-1./3.) - 25. * std::pow(double(A),-2./3.);
        
        //        double b=1.01*1.01*pow(double(A),1./3.);//OscillatorLength squared
        double b = phy::HBARC* phy::HBARC / ( hw * 938.9186945); //mass defined as avverage of proton neutron.
        return pow(b,L)*CoSMoValue;//convert to (e fm)^2L and divide by Weisskopf Estimate.
    }
    else if (type=='M')
    {
        double hw = 45. * std::pow(double(A),-1./3.) - 25. * std::pow(double(A),-2./3.);
        //        double b=1.01*1.01*pow(double(A),1./3.);//OscillatorLength squared
        double b = phy::HBARC* phy::HBARC /( hw * 938.9186945); //mass defined as avverage of proton neutron.
        return pow(b,L-1)*CoSMoValue;
    }
    //If this point is reached something went horribly wrong.
    cerr<<"Something is wrong in Weisskopf routine"<<endl;
    return -1;
}

double OscillatorLength(int A) {
    double hw = 45. * std::pow(double(A),-1./3.) - 25. * std::pow(double(A),-2./3.);
    
    //        double b=1.01*1.01*pow(double(A),1./3.);//OscillatorLength squared
    double b = phy::HBARC* phy::HBARC / ( hw * 938.9186945); //mass defined as avverage of proton neutron.
    return sqrt(b);
}

/*-----------MAIN---------------------------------------*/
int main (int argc, char** argv){
av::argparser args(argc,argv);
args.SetDescription("Analyze Electric/Magnetic/Beta transitions.");
args.SetDefault("XML file containing a parent system description.",2,2);
args.SetFlag("f","final","XML file containg a final daughter system, use if final system is in a separate xml",0,2);
args.SetFlag("o","output","Output XML file, if provided the default input is unchanged.",0,2);
args.SetFlag("w","Weisskopf","set to print Weisskopf units for E/M transitions. ",0,1);
    args.SetFlag("p","phase","set to change created operator's phase convention from R(0+ε) >0 to R(inf)>0 where ε -> 0+. ",0,1);
    args.SetFlag("u","units","set to print B(λΜ) in standard units:\n\t[>] e^2 (fm^(2L))\t for λ=Ε\n\t[>] μΝ^2 (fm)^(2L-2)\t for λ=Μ ",0,1);

    args.SetFlag("ft","ft"," print ft value for GT transition ",0,1);


if(args.ParseArguments()) return 1;
  //unsigned char dfd;
  /*----------Check that interaction file is fine----------*/
  string sysfile=args(1);
  XMLNode xMainNode=XMLNode::openFileHelper(sysfile.c_str(), xsm::cosmo); //we start by opening xml input file
  XMLNode xValence=xMainNode.getChildNode(xsm::valencespace);
  int levels=xValence.nChildNode(xsm::valencespace_orbital);
  //cerr<<levels<<endl;
 SM::SingleParticleStates sps; Make(sps,xValence);

  int systems=xMainNode.nChildNode(xsm::system);
//for testing purposes set first system
  //XMLNode xParentSys=xMainNode.getChildNode(xsm::system,-1);
  cerr<<"Select parent system \n";
  XMLNode xParentSys=xml::SelectNode(xMainNode,xsm::system,"name");
  
  XMLNode xMainNodeDaughter=xMainNode; // daughter system
  if (args.isSet("f",1)) {
  xMainNodeDaughter=XMLNode::openFileHelper(args("f",1).c_str(), xsm::cosmo);
  }
  cerr<<"Select daughter system \n";
  XMLNode xDaugterSys=xml::SelectNode(xMainNodeDaughter,xsm::system,"name");
  //XMLNode xDaugterSys=xMainNode.getChildNode(xsm::system,-1);
//Next we add operator

 
 SM::Many_Body_States stp;
 if(!stp.Read(xParentSys.getChildNode(xsm::system_basisstates).getAttribute("name"))) {FatalError("Many-body states (*.mbs) must be avalable");}
 SM::Many_Body_States std;
 if(!std.Read(xDaugterSys.getChildNode(xsm::system_basisstates).getAttribute("name"))) {FatalError("Many-body states (*.mbs) must be avalable");}

 //int levels=xValence.nChildNode(xsm::valencespace_orbital);
 

 
 //make set of quantum numbers for parent/daughter
 SM::State PT, DT;
 ReadXMLQN(PT, xParentSys);
 ReadXMLQN(DT, xDaugterSys);
 int kappa=-(PT.Jz-DT.Jz)/2;
 
 //operators
    bool phaseConvention = true;
    if (args.isSet("p"))
        phaseConvention = false;
    
 XMLNode xOperator=InquireXMLOperator(xMainNode); //end making or selecting operator
 MBOperator MUME(2);

 int K=ReadXMLOperator(MUME, kappa, sps, xOperator,phaseConvention);
 char tname[256];
 getXMLData(tname,"name",xOperator);
 
 

 
 //ElectricSPMultipole(MUME, sps, K, kappa, Rfi, charge[1], charge[0]);
 SM::PPMatrix ME(std,stp);
 ME.push_back(&MUME); 
 //MBOperator MUMM(2);
 //MagneticSPMultipole(MUMM, sps, K, kappa, Rfi);
 
 //MBOperator MUMGT(2);
 //BetaJAMultipole(MUMGT,sps,K,kappa,Rfi,0);
 //end operators
 
 
 //------------->operator for sum rule (square constructed directly<-----------
// tav::STensor<spsint,double> SPQ2(4); //two-body
// tav::STensor<spsint,double> SPQ1(2); //one-body of square
// cout<<K<<endl;
// //return 1;
//   for (int kappa=-K;kappa<=K;kappa++) {
//  JO.clear();
//  ReadXMLOperator(JO, kappa, sps, xOperator); 
////SM::StripHermitianOperator(JO); //many operators contain an extra conjugated part A+A^\dagger which never works
//  SM::SPSqrM(SPQ1,SPQ2,JO); //Square this operator
//  }
//  if(JO.empty()) {CriticalStop("Operator is empty");};
//  SM::StripHermitianOperator(SPQ1);
// SM::StripHermitianOperator(SPQ2);
// SM::PPMatrix MESQR(stp);
// MESQR.push_back(&SPQ1);
// MESQR.push_back(&SPQ2);
//  tav::Tensor<1,double> z(stp.n); //temporary vector
 //---------------><------------------------------
 
//prepare to go over all states in parent-daugther systems
 int pstates=xParentSys.nChildNode(xsm::system_eigenstate);
if (pstates==0) FatalError("no states in xml");
int dstates=xDaugterSys.nChildNode(xsm::system_eigenstate);
if (dstates==0) FatalError("no states in xml");
 XMLNode xParentStateNode;
 XMLNode xDaugterStateNode;
//    int Anucleus=1;
int Anucleus=1;
try {
getXMLData(Anucleus,"A",xParentSys);
//Anucleus=atoi(xParentSys.getAttribute("A"));
}
catch (...) {cerr<<"Nuclear mass number A is not defined, units will not work "; Anucleus=1;}
 ReadXMLQN(PT, xParentSys); //read system QN
 ReadXMLQN(DT, xDaugterSys); //read system QN
 
 for (int p=0;p<pstates;p++) {
 
 double bsum=0.0;
 cout<<endl;
 xParentStateNode=xParentSys.getChildNode(xsm::system_eigenstate,p);
  // cout<<xParentStateNode.createXMLString(true);
   ReadXMLState(PT, xParentStateNode); //read the rest of qn
//--------------------> second method <---------------------
//   MESQR.TimesVector(z,PT.z,0.0);
//   double xxsum=z*PT.z; //sum rule
//------------------>Alternative way of SUM RULE<------------------------------
//    double xxsum=0.0;
//for (int kappa=-K;kappa<=K;kappa++) {
// MBOperator JO(2);
//  ReadXMLOperator(JO, kappa, sps, xOperator); 
//SM::PPMatrix ME2(stp);
//ME2.push_back(&JO);
// xxsum+=ME2.Sqr(PT.z);
//  }
//----------------><-------------------------
   
 //----------------->Add sumrule to parent xml<---------------------------------
//     XMLNode xSumRule=xParentStateNode.addChild("sumrule");
//   char qname[20]; sprintf(qname,"B(%s)",tname);
//          addAttribute(xSumRule,"name", qname);
//	   addAttribute(xSumRule,"B",xxsum);
	   
//------------------>end add sumrule to parent xml<-----------------------
     for (int q=0;q<dstates;q++) {
		xDaugterStateNode=xDaugterSys.getChildNode(xsm::system_eigenstate,q);
	//	cout<<xDaugterStateNode.createXMLString(true);
		ReadXMLState(DT, xDaugterStateNode); //read the rest of qn
		string DaugterName; 
		getXMLData(DaugterName,"name",xDaugterStateNode);
        string DaugterE; 
		getXMLData(DaugterE,"E",xDaugterStateNode);
		//ReadXMLAttribute(DaugterName,"name");
   //  cout<<"-->xx "<<DT.E<<" "<<DT.J<<" ";
     int Jz=PT.Jz-DT.Jz;
     int Tz=PT.Tz-DT.Tz;
     //double htmp=1.0;
 
      double htmp=tav::ClebschGordan( DT.J, DT.Jz, 2*K, Jz, PT.J, PT.Jz);
       if (tav::Abs(htmp)<1E-10) continue;
       double xx=ME.OverlapNH(DT.z,PT.z);
       double BB=tav::Sqr(xx/htmp); //value of B(the reduced transition rate)
      // cerr<<tav::TThreeJSymbol(DT.J, DT.Jz, 2*K, -Jz, PT.J, PT.Jz)<<endl;
       //cerr<<DT.J<<" "<<DT.Jz<<" "<<2*K<<" "<< Jz<<" "<< PT.J<<" "<<PT.Jz<<endl;
       double MM=xx/tav::TThreeJSymbol(DT.J, DT.Jz, 2*K, -Jz, PT.J, PT.Jz); // we use our three j symbol with time reversed first argument (final state)
         double ww=-1,beu = -1;
       bsum+=BB;
           cout<<DT.J<<" "<< DT.Jz<<" "<< 2*K<<" "<< Jz<<" "<<PT.J<<" "<<PT.Jz<<"  M="<<MM<<" B="<<BB<<"  sum="<<bsum<<flush;
         if (args.isSet("w"))
         {
             ww=WeisskopfEstimate(xOperator,Anucleus,BB);
             if (ww>-1)
                 cout<<" W="<<ww;
         }
         if (args.isSet("u"))
         {
             beu = BLUnitConversion(xOperator,Anucleus,BB);
             if (beu > -1)
                 std::cout << " Bu=" << beu;
         }
        cout<<endl;
        
             
	   // At this stage we can compute transition
	   /*
       cout<<"--> "<<DT.E<<" "<<DT.J/2;
       cout<<" B("<<tname<<K<<")="<<tav::Sqr(xx/htmp);
       if (tname=='G') cout<<"  (ft)= "<<(6250./1.5129/(tav::Sqr(xx/htmp)))<<" ";
	   */
	   
//here is output for transition
	   XMLNode xTransition=xParentStateNode.addChild("transition");
	   char qname[20]; sprintf(qname,"B(%s)",tname);
	  // cerr<<qname<<endl;
       addAttribute(xTransition,"name", qname);
	  // addAttribute(xTransition,xsm::J,spin(DT.J));
	  // addAttribute(xTransition,xsm::T,spin(DT.T));
	    addAttribute(xTransition,"final",DaugterName.c_str());
         addAttribute(xTransition,"E",DaugterE.c_str());
		if (tname[0]=='G') {
		addAttribute(xTransition,"units","ga^2/(4pi)");
	  }
	   addAttribute(xTransition,"B",BB);
	   if (args.isSet("w") && ww>-1) addAttribute(xTransition,"weisskopf",ww);
       if (args.isSet("u") && beu>-1) addAttribute(xTransition,"Bu",beu);
      // addAttribute(xTransition,"fraction",BB/xxsum);
	   //see BMI, page 349, Suhonen 172
     //ZV 24.33
         
         if (args.isSet("ft")) {
             double TT=6145.0;
             double La=1.27;
addAttribute(xTransition,"logft",log(TT/La/La/(BB))/log(10.0));
         }
        
        //output MOMENT-Diagonal transition
		 if ((fabs(PT.E-DT.E)<1E-5)&&(PT.J==DT.J)&&(PT.T==DT.T)&&(K<3)&&((PT.Tz==DT.Tz)))   {
		 XMLNode xTransition=xParentStateNode.addChild("moment", false, 0);
		 addAttribute(xTransition,"name", tname);
         double mmx=xx/htmp*tav::ClebschGordan(DT.J, DT.J, 2*K, 0, PT.J, PT.J)*sqrt(4.0*PI/(2.0*K+1))*K;
		 addAttribute(xTransition,"m",mmx);
         if (args.isSet("u")) addAttribute(xTransition,"oscl",OscillatorLength(Anucleus));
		 continue;
		 }
		 
		 
        //note no 4\pi because my operator has no 4\pi
       //if ((fabs(PT.E-DT.E)<1E-5)&&(PT.J==DT.J)&&(PT.T==DT.T)) {
	    //cout<<" <M>= "<<xx/htmp*tav::ClebschGordan(DT.J, DT.J, 2*K, 0, PT.J, PT.J);
	//cout<<" moment(E1,E2,M1,M2)= "<<xx/htmp*tav::ClebschGordan(DT.J, DT.J, 2*K, 0, PT.J, PT.J)*sqrt(4.0*PI/(2.0*K+1))*K;
	   //note we multiply by K so it works for K=1 and K=2 as deformation moment  
       }
       //cout<<" Q="<<PT.E-DT.E<<endl;
       //add fraction of sumrule in transitions computed
	   //addAttribute(xSumRule,"fraction",bsum/xxsum);
   }
   

 
 //cout<<xMainNode.createXMLString(true)<<endl;
 //if (argc==3) xMainNode.writeToFile(argv[2]);

if (args.isSet("o",1)) xMainNode.writeToFile(args("o",1).c_str());
else xMainNode.writeToFile(sysfile.c_str());
 return 0;
}

  







