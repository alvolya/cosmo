/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/** 
\author Alexander Volya  <http://www.volya.net>
\file sphericalint.cxx
\brief Control and Manipulation with two-body interactions in spherical basis

\todo Change System to Space.
*/
/*
06/30/09 Drop all dependence on System/Space
01/18/10 ReadOxbashINT reads SMINTPATH variable, bug fix with comma in s.p.e
*/

//#include "System.cxx"
#ifndef __SPHERICALINT_CXX_
#define __SPHERICALINT_CXX_
#include <cstdlib> 
using std::getenv;
#include <inquire.cxx>
using tav::Inquire;
//use getenv to read path of the interaction files
#include <av/iofile.cxx> //reading file with path
#include <tav/basic.cxx>
using tav::Min;
using tav::Max;
using tav::Parity;
#include <tav/stensor.cxx>
#include <tav/tensor.cxx>
#include <mav/SixJSymbol.cxx>
#include <av/TextAnalysis.cxx> //Needed for ReadSphericalINT
#include <regex> //line counts
             //NextFermiOperators
using mav::SixJSymbol;
namespace SM {
 /*! \brief Check spherical int for hermicity and cut dublicate parts
 */
   int CheckSphericalINT(tav::STensor<int,double> &INT, bool zerobad=true) {
   int ret0=0;
   tav::STensor<int,double>::iterator ip1,ip2;
    for (ip1=INT.begin();ip1!=INT.end();ip1++) {
		//read data in sorted order by pair
	int a[4]; 
	if ((ip1->first[0])<(ip1->first[1])) {a[0]=(ip1->first)[0]; a[1]=(ip1->first)[1];} 
	else {a[1]=(ip1->first)[0]; a[0]=(ip1->first)[1];}
	if ((ip1->first[2])<(ip1->first[3])) {a[2]=(ip1->first)[2]; a[3]=(ip1->first)[3];} 
	else {a[2]=(ip1->first)[3]; a[3]=(ip1->first)[2];}
	ip2=ip1;
	++ip2;
	for (;ip2!=INT.end();ip2++)
	{
	if (ip1->first[4]!=ip2->first[4]) continue; //spins must be the same
	if (ip1->first[5]!=ip2->first[5]) continue; //isospins must be the same
	int b[4]; 	
	if ((ip2->first[0])<(ip2->first[1])) {b[0]=(ip2->first)[0]; b[1]=(ip2->first)[1];} 
	else {b[1]=(ip2->first)[0]; b[0]=(ip2->first)[1];}
	if ((ip2->first[2])<(ip2->first[3])) {b[2]=(ip2->first)[2]; b[3]=(ip2->first)[3];} 
	else {b[2]=(ip2->first)[3]; b[3]=(ip2->first)[2];}
	bool badline;
	badline=((a[0]==b[0])&&(a[1]==b[1])&&(a[2]==b[2])&&(a[3]==b[3]));
	badline=badline||((a[0]==b[2])&&(a[1]==b[3])&&(a[2]==b[0])&&(a[3]==b[1]));
	if (badline) {
//	cerr<<"WARNING: There is repeated interaction (your may double count!!!)  "<<endl;
//	cerr<<ip1->first[0]<<" "<<ip1->first[1]<<" "<<ip1->first[2]<<" "<<ip1->first[3]<<" "<<ip1->first[4]<<" "<<ip1->first[5]<<" "<<ip1->second<<endl;
//	cerr<<ip2->first[0]<<" "<<ip2->first[1]<<" "<<ip2->first[2]<<" "<<ip2->first[3]<<" "<<ip2->first[4]<<" "<<ip2->first[5]<<" "<<ip2->second;
//	zerobad=true; //do not ask questions
if (ret0==0) {cerr<<"There are repeated components "<<endl;
if (zerobad) cerr<<"These components are removed "<<endl;
}
//    if (!zerobad) {cerr<<endl; zerobad=Inquire("Zero this and all following repeted components ", true);}
	if (zerobad) {(ip2->second)=0.0;}
    ret0++;
	//INT.erase(ip2);//we assume that now iterator points to the next element
	}
	//else ++ip2; //increment element
	} //loop ip2
   } //loop ip1
   return ret0;
   };
   
  /*! \brief read Interaction in oxbash format from stream
  ignore single-particle energies by skipping line 
  ignore all scaling defined in oxbash format
   */
  int ReadSphericalINT(
      tav::STensor<int,double> &INT, ///<output: sparse-tensor with data V[1][2][3][4][L][T] 
  std::ifstream &DataStream ///<input: stream
                      ) 
  {
    const int MaxLineLength=1024;
    char   strFileLine[MaxLineLength]; // One line buffer...
    long FPosition;
//    DataStream.seekg(std::ios::beg); //this line is important
    FPosition=av::SkipComments(DataStream);
    if (FPosition==-1) FatalError("No data is available in stream");    
    int intlines=0;
    DataStream>>intlines;  //in testing version skip the rest
    DataStream.getline(strFileLine, MaxLineLength, '\n');
    intlines=abs(intlines);
    int ijk[6]; 
    double htmp;
    for (int il=0;il<intlines;il++) {
      DataStream>>ijk[2];
      DataStream>>ijk[3];
      DataStream>>ijk[0];
      DataStream>>ijk[1];
      DataStream>>ijk[4];
      DataStream>>ijk[5];
      DataStream>>htmp;
	  CheckInputStream(DataStream);
  //DataStream>>VX;
  //convert to from oxbash notation
      ijk[0]--;
      ijk[1]--;
      ijk[2]--;
      ijk[3]--;
  // Note OXBASH ORDERING IS FROM HIGHER TO LOWER 
 // if (ijk[0]>ijk[1]) tav::Swap(ijk[0],ijk[1]);
 // if (ijk[2]>ijk[3]) tav::Swap(ijk[2],ijk[3]);
 // if (ijk[0]>ijk[2]) {tav::Swap(ijk[0],ijk[2]); tav::Swap(ijk[1],ijk[1]);
      INT[ijk]+=htmp;
    }
    return intlines;
  };
/** \brief Read Oxbash-format interaction from file
 */ 
  int ReadSphericalINT(tav::STensor<int,double> &INT, const char *sysfile)
  {

    std::ifstream DataStream;  
    DataStream.open(sysfile);
    if (!DataStream)   {               // if the file does not exist
      FatalError("ReadSphericalINT: File "<<sysfile<<" is not found \n"); }
      ReadSphericalINT(INT, DataStream);
      DataStream.close();
      return 1;
  }
   /*! \brief read OXBASH (*.int) interaction from stream 
   */
  int ReadOxbashINT(
      std::vector<double> &ENT,
     // tav::STensor<int,double> &ENT, ///<output: vector data e[1][2]
      tav::STensor<int,double> &INT, ///<output: sparse-tensor with data V[1][2][3][4][L][T] 
	  std::ifstream &DataStream, ///<input: stream
	  int N, ///number of particles for scaling
	  double spescaling, ///<additional scaling of s.p. energies
	  double vscaling ///<additional scaling of matrix elements (besides those in file)
                      ) 
  {   //prepare data
      int ijk[6];
       tav::STensor<int,double> INTTMP(6); //two-body
      //read file line-by line and run regex
      int linecount=0; //count non-comment lines
      std::string FileLine;
      int intlines=0; //lines as read from file
      while (std::getline(DataStream, FileLine)) {
              if (FileLine.empty()) continue; //check for empty string
              //std::regex re_empty("^\\s*");  //any number of spaces \\s* followd by 0 or 1  of the special characters
             //using regex to detect comment line does not work ??? use old method
             //->> if(std::regex_match(FileLine, std::regex("^\\s*[!#$%]?"))) continue; // string has only spaces
              FileLine=std::regex_replace(FileLine, std::regex("^\\s*"),""); //replace head spaces with nothing
              if (av::CommentString(FileLine.c_str())) continue;
              //------use regex to check comment -does not seem to fully work---
             // if(std::regex_match(FileLine, std::regex("\\s*[!#%].*"))) continue; // string has only spaces
             // if (FileLine.empty()) continue; //check again if empty
              //------end regex to check comment
              //std::cout<<FileLine<<std::endl;
              
              //std::regex re("[\\|,:]");
              //std::regex re("\\s+");
              // Delimiters are spaces (\s) and/or commas
              std::regex re("[\\s,]+");
              std::sregex_token_iterator first{FileLine.begin(), FileLine.end(), re, -1}, last;//the '-1' is what makes the regex split (-1 := what was not matched)
              std::vector<std::string> tokens{first, last};
//              std::cout<<"Line number :"<<linecount<<endl;
//              for (int i=0;i<tokens.size();++i) std::cout<<tokens[i]<<" ";
//              std::cout<<tokens.size()<<endl;
			//--------------Here we have our line-----------------------
          if (linecount==0) {
              //we are dealing with the first line
              int words=tokens.size(); //number of words in the first line
              //if (is_int(tokens[0]))
                  intlines=std::stoi(tokens[0]);
              //else {cerr<<"Error reading number of lines\n"; return 1;}
              words--; //discount number of lines
              if (intlines<0) words-=3; //last 3 are scaling
              if (intlines==0) cerr<<"The number of interactions to read should not be zero ?"<<endl;
              //CriticalError("The number of interactions to read should not be zero");
              
             // double htmp;
              //
              if (ENT.size()<words) ENT.resize(words); //resize elements
              //
              for (int i=0;i<words;i++) {
              //if (is_float(tokens[i+1]))
                  ENT[i]+=std::stof(tokens[i+1])*spescaling;
                  //else {cerr<<"Error reading single particle energy \n"; return 1;}
              }
              //DataStream.getline(strFileLine, MaxLineLength, '\n');
                double ACORE; //scaling=(ASCALE/(N+ACORE)^SPOW
                double ASCALE;
                double SPOW;
              if (intlines<0) {
                  ACORE=std::stof(tokens[words+1]);
                  ASCALE=std::stof(tokens[words+2]);
                  SPOW=std::stof(tokens[words+3]);
                if (N>=0) {vscaling*=std::pow(ASCALE/(N+ACORE),SPOW);//we add this scaling on top
                //cerr<<"We have some scaling "<<vscaling<<endl;
                }
              }
              intlines=abs(intlines);
              INF(cerr<<"Single-particle levels    :" << words <<" scale: "<<spescaling<<" \n");
              INF(cerr<<"Lines in interaction file :" << intlines <<" scale: "<<vscaling<<" \n");
              
          } //end first line read
          else {
              if (tokens.size()<7) continue; //make sure there is at least 7 entries in line
              //for (int i=0;i<6;++i) cerr<<tokens[i]<<" ";
              // cerr<<endl;
                   ijk[2]=std::stoi(tokens[0]);
                   ijk[3]=std::stoi(tokens[1]);
                   ijk[0]=std::stoi(tokens[2]);
                   ijk[1]=std::stoi(tokens[3]);
                   ijk[4]=std::stoi(tokens[4]);
                   ijk[5]=std::stoi(tokens[5]);
 
              //convert to from oxbash notation
              		ijk[0]--;
             	 	ijk[1]--;
              		ijk[2]--;
              		ijk[3]--;
                   INTTMP[ijk]+=vscaling*std::stof(tokens[6]);
          }//end read other lines
          linecount++;
          }//end loop over file
      if ((intlines!=0)&&(abs(intlines)!=linecount-1)) cerr<<"ReadOxbashINT: Line number mismatch read "<<linecount-1<<" listed in file "<<intlines<<endl;

	CheckSphericalINT(INTTMP); //conduct check
	INT+=INTTMP; //use add tensor method
    return intlines;
  };
//this is old string-based reader requires interaction file not to have extra terms
 int ReadOxbashINT1(
      std::vector<double> &ENT,
     // tav::STensor<int,double> &ENT, ///<output: vector data e[1][2]
      tav::STensor<int,double> &INT, ///<output: sparse-tensor with data V[1][2][3][4][L][T]
      std::ifstream &DataStream, ///<input: stream
      int N, ///number of particles for scaling
      double spescaling, ///<additional scaling of s.p. energies
      double vscaling ///<additional scaling of matrix elements (besides those in file)
                      )
  {
    const int MaxLineLength=1024;
    char   strFileLine[MaxLineLength]; // One line buffer...
    long FPosition;
//    DataStream.seekg(std::ios::beg); //this line is important
    FPosition=av::SkipComments(DataStream);
    if (FPosition==-1) FatalError("No data is available in stream");
//just count lines and read them
    FPosition=DataStream.tellg();
    DataStream.getline(strFileLine, MaxLineLength, '\n'); //read next line
    int words=av::WordCount(strFileLine);
    DataStream.seekg(FPosition); //return back
//end count lines
      
    int intlines=0;
    DataStream>>intlines;  //in testing version skip the rest
    words--;
    if (intlines<0) words-=3; //last 3 are scaling
    if (intlines==0) cerr<<"The number of interactions to read should not be zero ?"<<endl;
    //CriticalError("The number of interactions to read should not be zero");
    
    double htmp;
    //
    if (ENT.size()<words) ENT.resize(words); //resize elements
    //
    for (int i=0;i<words;i++) {
    DataStream>>htmp; ENT[i]+=htmp*spescaling;
    if (DataStream.peek()==',') DataStream.ignore(); //ignore "," in some some oxbash data files
    }
    //DataStream.getline(strFileLine, MaxLineLength, '\n');
      double ACORE; //scaling=(ASCALE/(N+ACORE)^SPOW
      double ASCALE;
      double SPOW;
    if (intlines<0) {
      DataStream>>ACORE;
      DataStream>>ASCALE;
      DataStream>>SPOW;
      if (N>=0) {vscaling*=std::pow(ASCALE/(N+ACORE),SPOW);//we add this scaling on top
      //cerr<<"We have some scaling "<<vscaling<<endl;
      }
    }
    intlines=abs(intlines);
    INF(cerr<<"Single-particle levels    :" << words <<" scale: "<<spescaling<<" \n");
    INF(cerr<<"Lines in interaction file :" << intlines <<" scale: "<<vscaling<<" \n");
    int ijk[6];
    //read first into temporary file
     tav::STensor<int,double> INTTMP(6); //two-body
    for (int il=0;il<intlines;il++) {
      DataStream>>ijk[2];
      DataStream>>ijk[3];
      DataStream>>ijk[0];
      DataStream>>ijk[1];
      DataStream>>ijk[4];
      DataStream>>ijk[5];
      DataStream>>htmp;
       CheckInputStream(DataStream);
      DBG(cerr<<ijk[0]<<" "<<ijk[1]<<" "<<ijk[2]<<" "<<ijk[3]<<" "<<ijk[4]<<" "<<ijk[5]<<" "<<htmp<<endl);
  //DataStream>>VX;
  //convert to from oxbash notation
      ijk[0]--;
      ijk[1]--;
      ijk[2]--;
      ijk[3]--;
  // Note OXBASH ORDERING IS FROM HIGHER TO LOWER
 // if (ijk[0]>ijk[1]) tav::Swap(ijk[0],ijk[1]);
 // if (ijk[2]>ijk[3]) tav::Swap(ijk[2],ijk[3]);
 // if (ijk[0]>ijk[2]) {tav::Swap(ijk[0],ijk[2]); tav::Swap(ijk[1],ijk[1]);
      INTTMP[ijk]+=htmp*vscaling;
    }
    CheckSphericalINT(INTTMP); //conduct check
    INT+=INTTMP; //use add tensor method
    return intlines;
  };
  //compatibility 
    int ReadOxbashINT(
      tav::STensor<int,double> &ENT, ///<output: vector data e[1][2]
      tav::STensor<int,double> &INT, ///<output: sparse-tensor with data V[1][2][3][4][L][T] 
	  std::ifstream &DataStream, ///<input: stream
	  int N, ///number of particles for scaling
	  double spescaling, ///<additional scaling of s.p. energies
	  double vscaling ///<additional scaling of matrix elements (besides those in file)
	)
	{
	std::vector<double> EE;
	int intlines=ReadOxbashINT(EE, INT, DataStream, N, spescaling, vscaling );
	for (int i=0;i<EE.size();i++) ENT(i,i)+=EE[i];
	return intlines;
	} 
/** \brief Read Oxbash-format interaction from file
 */ 
   template <typename myspvector>
  int ReadOxbashINT(
      myspvector &ENT, ///<output: vector data e[1]
      tav::STensor<int,double> &INT, ///<output: sparse-tensor with data V[1][2][3][4][L][T] 
	  const char *sysfile,
	  int N,
	  double spescaling=1.0, ///<additional scaling of s.p. energies
	  double vscaling=1.0 ///<additional scaling of matrix elements (besides those in file)
	  )
  {
    std::ifstream DataStream;  
    DataStream.open((av::FileName(sysfile,"SMINTPATH",".INT")).c_str());
      //if (!DataStream)   {               // if the file does not exist
      //FatalError("ReadSphericalINT: File "<<sysfile<<" is not found \n"); }
	  INF(std::cerr<<"Reading Interactions from \n"<<av::FileName(sysfile,"SMINTPATH",".INT")<<std::endl);
      ReadOxbashINT(ENT,INT, DataStream,N,spescaling,vscaling);
      DataStream.close();
      return 1;
  }

  /**!\brief Generate a list of all possible two-body interactions
  */
  int ListINT(ValenceSpace &x) {
  int s1,s2,s3,s4;

    if (x.size()==0) FatalError("ListINT: you must define the system first");
	int levels=x.size();
  //int sg=1-2*x.type; //-1 fermions +1 bosons
  //we first need to count how many pair states we can have 
  //cerr<<levels<<endl;
  int tbme=0;
  for (s2=0;s2<levels;s2++) 
    for (s1=s2;s1<levels;s1++) //this will list full pair (s1>= s2) -oxbash  
      for (s4=s2;s4<levels;s4++) {
	if (s4==s2) s3=Max(s4,s1); else s3=s4;
	//this is where the trick is we compare pair by first number s2 s4 
	//but if those are equal we compare using s1 and s3
	//rmemember s1>=s2 same s3>=s4 
	//from independence of pairs 
	//s4>=s2 and if s4==s2 then s3>=s1
	for (;s3<levels;s3++){
	  //  cerr<<"("<<s1<<" "<<s2<<" )  ("<<s3<<" "<<s4<<" )";
       //now we have all pairs and we need to decide what L and T are alowed
	    for (int L=abs((x[s1].J-x[s2].J)/2);L<=((x[s1].J+x[s2].J)/2);L++) 
	      for (int T=abs((x[s1].T-x[s2].T)/2);T<=((x[s1].T+x[s2].T)/2);T++)
  {
	      //now we need to check if L and t are alowed for both pairs
	      if ((2*L)>(x[s3].J+x[s4].J)) continue; //bad triangle rule
	      if ((2*L)<abs(x[s3].J-x[s4].J)) continue; //bad triangle rule
	      if ((2*T)>(x[s3].T+x[s4].T)) continue; //bad triangle rule
	      if ((2*T)<abs(x[s3].T-x[s4].T)) continue; //bad triangle rule
	      //now if we are on one level I need to check parity
	      if ((s1==s2)&&(((x[s1].J-L+x[s1].T-T)%2)==0)) continue;
              if ((s3==s4)&&(((x[s3].J-L+x[s3].T-T)%2)==0)) continue;
	       if (((x[s1].L+x[s2].L+x[s3].L+x[s4].L)&1)) continue; //parity
	      tbme++; 
	      cerr<<"("<<s1<<" "<<s2<<" )  ("<<s3<<" "<<s4<<" )"<<" "<<L<<" "<<T<<endl;
  }//this is end cycle over T & L
	}}
//this goes over pair (s3>=s4) 
	    /*let us define pair comparison so not to double count pairs
	      (s1,s2) (s3,s4)  if s2>s4 then P12>P34
	
     */
  cerr<<tbme<<endl;
  DBGPrint(tbme);
  return tbme;
  }
/** \brief Generate empty sparse tensor structure with restrictions set by 
valence space should work for fermions and bosons line 
      if ((s1==s2)&&(((x[s1].J-L+x[s1].T-T)%2)==(x[s1].J+1)%2)) continue;
controls what to skip
*/
  int GenerateEmptyINT(tav::STensor<int,double> &VSP, ValenceSpace &x) {
  int s1,s2,s3,s4;
  int *aa=new int [6];
 if (x.size()==0) FatalError("GenerateEmptyINT: you must define the system first");
	int levels=x.size();
	//cerr<<"we have levels "<<levels<<endl;
 // int sg=1-2*x.type; //-1 fermions +1 bosons
  //we first need to count how many pair states we can have 
 // cerr<<levels<<endl;
  int tbme=0;
  for (s2=0;s2<levels;s2++) 
    for (s1=s2;s1<levels;s1++) //this will list full pair (s1>= s2) -oxbash  
      for (s4=s2;s4<levels;s4++) {
	if (s4==s2) s3=Max(s4,s1); else s3=s4;
	//this is where the trick is we compare pair by first number s2 s4 
	//but if those are equal we compare using s1 and s3
	//rmemember s1>=s2 same s3>=s4 
	//from independence of pairs 
	//s4>=s2 and if s4==s2 then s3>=s1
	for (;s3<levels;s3++){
	  //  cerr<<"("<<s1<<" "<<s2<<" )  ("<<s3<<" "<<s4<<" )";
       //now we have all pairs and we need to decide what L and T are alowed
	    for (int L=abs((x[s1].J-x[s2].J)/2);L<=((x[s1].J+x[s2].J)/2);L++) 
	      for (int T=abs((x[s1].T-x[s2].T)/2);T<=((x[s1].T+x[s2].T)/2);T++)
  {
	      //now we need to check if L and t are alowed for both pairs
	      if ((2*L)>(x[s3].J+x[s4].J)) continue; //bad triangle rule
	      if ((2*L)<abs(x[s3].J-x[s4].J)) continue; //bad triangle rule
	      if ((2*T)>(x[s3].T+x[s4].T)) continue; //bad triangle rule
	      if ((2*T)<abs(x[s3].T-x[s4].T)) continue; //bad triangle rule
	      //now if we are on one level I need to check parity
          //(x[s1].J+1)%2 gives 0 for fermions 1 for bosons
	      if ((s1==s2)&&(((x[s1].J-L+x[s1].T-T)%2)==(x[s1].J+1)%2)) continue;
          if ((s3==s4)&&(((x[s3].J-L+x[s3].T-T)%2)==(x[s2].J+1)%2)) continue;
	      if (((x[s1].L+x[s2].L+x[s3].L+x[s4].L)&1)) continue; //parity
	      tbme++; 
	      // cerr<<"("<<s1<<" "<<s2<<" )  ("<<s3<<" "<<s4<<" )"<<" "<<L<<" "<<T<<endl;

	      aa[0]=s1; aa[1]=s2; aa[2]=s3; aa[3]=s4; aa[4]=L; aa[5]=T;
	      VSP[aa]=0.0;

  }//this is end cycle over T & L
	}}
//this goes over pair (s3>=s4) 
	    /*let us define pair comparison so not to double count pairs
	      (s1,s2) (s3,s4)  if s2>s4 then P12>P34
	
     */
  // cerr<<tbme<<endl;
  // DBGPrint(tbme);
  delete [] aa;
  return tbme;
  }
/** \brief Generate empty sparse tensor for quadrupole-quadrupole interaction
*/
 int GenerateEmptyQINT(tav::STensor<int,double> &VSP, ValenceSpace &x) {
  int s1,s2,s3,s4;
  int *aa=new int [6];

    if (x.size()==0) FatalError("GenerateEmtpyQINT: you must define the system first");
	int levels=x.size(); //we first need to count how many pair states we can have 
  //cerr<<levels<<endl;
  int tbme=0;
  for (s2=0;s2<levels;s2++) 
    for (s1=s2;s1<levels;s1++) //this will list full pair (s1>= s2) -oxbash  
      for (s4=s2;s4<levels;s4++) {
	if (s4==s2) s3=Max(s4,s1); else s3=s4;
	//this is where the trick is we compare pair by first number s2 s4 
	//but if those are equal we compare using s1 and s3
	//rmemember s1>=s2 same s3>=s4 
	//from independence of pairs 
	//s4>=s2 and if s4==s2 then s3>=s1
	for (;s3<levels;s3++){
	  //  cerr<<"("<<s1<<" "<<s2<<" )  ("<<s3<<" "<<s4<<" )";
       //now we have all pairs and we need to decide what L and T are alowed
	    for (int L=abs((x[s1].J-x[s2].J)/2);L<=((x[s1].J+x[s2].J)/2);L++) 
	      for (int T=abs((x[s1].T-x[s2].T)/2);T<=((x[s1].T+x[s2].T)/2);T++)
  {
	      //now we need to check if L and t are alowed for both pairs
	      if ((2*L)>(x[s3].J+x[s4].J)) continue; //bad triangle rule
	      if ((2*L)<abs(x[s3].J-x[s4].J)) continue; //bad triangle rule
	      if ((2*T)>(x[s3].T+x[s4].T)) continue; //bad triangle rule
	      if ((2*T)<abs(x[s3].T-x[s4].T)) continue; //bad triangle rule
	      //now if we are on one level I need to check parity
	      //if ((s1==s2)&&(((x[s1].J-L+x[s1].T-T)%2)==0)) continue;
              //if ((s3==s4)&&(((x[s3].J-L+x[s3].T-T)%2)==0)) continue;
	       if (((x[s1].L+x[s2].L+x[s3].L+x[s4].L)&1)) continue; //parity
	      tbme++; 
	      // cerr<<"("<<s1<<" "<<s2<<" )  ("<<s3<<" "<<s4<<" )"<<" "<<L<<" "<<T<<endl;

	      aa[0]=s1; aa[1]=s4; aa[2]=s3; aa[3]=s2; aa[4]=L; aa[5]=T;
	      VSP[aa]=0.0;

  }//this is end cycle over T & L
	}}
//this goes over pair (s3>=s4) 
	    /*let us define pair comparison so not to double count pairs
	      (s1,s2) (s3,s4)  if s2>s4 then P12>P34
	
     */
  // cerr<<tbme<<endl;
  // DBGPrint(tbme);
  delete [] aa;
  return tbme;
  }

/** \brief Pandeya transformation
This is a transformation PP to QQ (add formt), second application of this operation reverses the transformation. 
  \f[  {\tilde{V}}_{KS}^{({\bf 1}{\bf 4};{\bf 3}{\bf 2})}=\sum_{LT}(2L+1)(2T+1) \,\left\{ \begin{array}{ccc}
  1/2 & 1/2 & T\\
  1/2 & 1/2 & S\end{array}\right\} \left\{ \begin{array}{ccc}
  j_{1} & j_{2} & L\\
  j_{3} & j_{4} & K\end{array}\right\} \, V_{L\, T}^{({\bf 1}{\bf 2};{\bf 3}{\bf 4})}\,,
  \f]
\ref p_pandeya
*/
  int QQ2PP(tav::STensor<int,double> &VSP, ///<output 
            tav::STensor<int,double> &VSQ, ///<input
            ValenceSpace &x) {
    tav::STensor<int,double>::iterator ip,iq;
    int *ap;
    int *aq;
    double htmp;
    for (ip=VSP.begin();ip!=VSP.end();ip++) {
      (*ip).second=0.0; //zero first
       for (iq=VSQ.begin();iq!=VSQ.end();iq++) 
      {
	if ((*ip).first[0]!=(*iq).first[0]) continue;
	if ((*ip).first[1]!=(*iq).first[3]) continue;
	if ((*ip).first[2]!=(*iq).first[2]) continue;
	if ((*ip).first[3]!=(*iq).first[1]) continue;
	ap=(*ip).first; //put a pointer to the beguinning 
	aq=(*iq).first; //put a pointer to the beguinning 
	htmp=SixJSymbol(x[ap[2]].J/2., x[ap[1]].J/2., double(aq[4]) ,
		   x[ap[0]].J/2., x[ap[3]].J/2., double(ap[4]))*
	SixJSymbol(x[ap[2]].T/2., x[ap[1]].T/2., double(aq[5]) ,
		   x[ap[0]].T/2., x[ap[3]].T/2., double(ap[5]));
	//	cout<<htmp<<endl;
	if (fabs(htmp)>1E-8) 
	  (*ip).second+=(2.0*aq[4]+1.0)*(2.0*aq[5]+1.0)*htmp*((*iq).second);
	  
      }}
    return 1;
  }

  /** \brief Transform (add) interaction  in sparse tensor into system format 
  \todo Factor of 2 is removed
  */
  int INTS2T(tav::Tensor<6,double> &V, ValenceSpace &x, tav::STensor<int,double> &VSP) {  
    int s1,s2,s3,s4,l,t;
    double htmp;
    int sg=1-2*(x[0].J%2); //-1 fermions +1 bosons
	//prepare tensor 
	spin  Lmax=0;
	spin  Tmax=0;
    for (int ll=0;ll<x.size();ll++) {
	if (Lmax<(x[ll].J)) Lmax=x[ll].J;
	 if (Tmax<(x[ll].T)) Tmax=x[ll].T;
    }
      //These two are here because the casting was not working properly
      //and gave malloc errors when casting from SM::spin
      unsigned long long lmax =Lmax +1;
      unsigned long long tmax = Tmax+1;
   V.FullResize(x.size(),x.size(),x.size(),x.size(),lmax,tmax);
   V=0.0;//gives bug if not defined
	//end prepare tensor
	
	
    tav::STensor<int,double>::iterator ip;
    for (ip=VSP.begin();ip!=VSP.end();ip++) {
      s1=(*ip).first[0];
      s2=(*ip).first[1];
      s3=(*ip).first[2];
      s4=(*ip).first[3];
      l=(*ip).first[4];
      t=(*ip).first[5];
      htmp=(*ip).second;

  //change for permutations
      //htmp*=(sqrt((1.0+delta(s1,s2))*(1.0+delta(s3,s4)))/2.0);
	  if (s1!=s2) htmp/=sqrt(2.0);
	  if (s3!=s4) htmp/=sqrt(2.0);

  /*when we add there is a problem as we add more times and I need to
      accout for repeated indexes 
      Below is a fix for extra terms
  */

      if (s1==s2) htmp/=2.0; //statistics in state in
      if (s3==s4) htmp/=2.0; //statistics in state out
      if ((s1==s3)&&(s2==s4)) htmp/=2.; //T invariance
  

      V[s1][s2][s3][s4][l][t]+=htmp;
 //T-invariant
      V[s3][s4][s1][s2][l][t]+=htmp; 
      V[s2][s1][s3][s4][l][t]+=sg*htmp*
          Parity((x[s1].J+x[s2].J-2*l+x[s1].T+x[s2].T-2*t)/2);
      V[s3][s4][s2][s1][l][t]+=sg*htmp*
          Parity((x[s1].J+x[s2].J-2*l+x[s1].T+x[s2].T-2*t)/2);
      V[s1][s2][s4][s3][l][t]+=sg*htmp*
          Parity((x[s3].J+x[s4].J-2*l+x[s3].T+x[s4].T-2*t)/2);
      V[s4][s3][s1][s2][l][t]+=sg*htmp*
          Parity((x[s3].J+x[s4].J-2*l+x[s3].T+x[s4].T-2*t)/2);
      V[s4][s3][s2][s1][l][t]+=
          htmp*Parity((x[s3].J+x[s4].J-2*l+x[s3].T+x[s4].T-2*t)/2)*
          Parity((x[s1].J+x[s2].J-2*l+x[s1].T+x[s2].T-2*t)/2);
      V[s2][s1][s4][s3][l][t]+=
          htmp*Parity((x[s3].J+x[s4].J-2*l+x[s3].T+x[s4].T-2*t)/2)*
          Parity((x[s1].J+x[s2].J-2*l+x[s1].T+x[s2].T-2*t)/2);

    }

    return 1;
  }
  
  template <typename myspvector>
  int ReadSminINT(
      myspvector &ENT, ///<output: vector data e[1]
      tav::STensor<int,double> &INT, ///<output: sparse-tensor with data V[1][2][3][4][L][T] 
	  std::ifstream &DataStream, ///<input: stream
	  int N ///<particle number
                      ) 
{
   int files;
   /*
   int N;
   if (ReadTextData(N,"N",DataStream)) {
   cerr<<"Particle number is not defined, enter N= ";
   cin>>N;
   }
*/
   if (ReadTextData(files,"files",DataStream)) CriticalError("There are no interaction files to read");
   char fname[128];
   double spescaling,vscaling;
   for (int i=0;i<files;i++) {
   DataStream>>fname;
   DataStream>>spescaling;
   DataStream>>vscaling;
   INF(cerr<<"Reading int file "<<fname<<" "<<spescaling<<" "<<vscaling<<"\n"); 
   ReadOxbashINT(ENT,INT,fname,N,spescaling,vscaling);
   }
	return 0;
}
  template <typename myspvector>
 int ReadSminINT(
      myspvector &ENT, ///<output: vector data e[1]
      tav::STensor<int,double> &INT, ///<output: sparse-tensor with data V[1][2][3][4][L][T] 
	  const char * sysfile, ///<input: stream
	  int N ///<particle number
                      )
					    {

    std::ifstream DataStream; 
	 
  //  DataStream.open(sysfile);
	DataStream.open((av::FileName(sysfile,"SMINTPATH",".smin")).c_str());
    if (!DataStream)   {               // if the file does not exist
      FatalError("ReadSminINT: File "<<sysfile<<" is not found \n"); }
      ReadSminINT(ENT,INT, DataStream,N);
      DataStream.close();
      return 0;
  }
 /*! \brief Read single particle energies
 */
   int ReadSPE(
       std::vector<double> &x, ///< Spherical, system data, x.e is used. 
	  std::ifstream &DataStream ///<input: stream
                      ) 
{
  int levels;
   if (ReadTextData(levels,"type",DataStream)) FatalError("We need to know the number of s.p. levels");
    int *e=new int [levels];
	x.resize(levels);
	if (ReadTextData(e,levels,"spenergies",DataStream)) 
	{cerr<<"SPE are not included \n"; for (int i=0;i<levels;i++) x[i]=0.0;}
	else {for (int i=0;i<levels;i++) x[i]=e[i];}
    delete [] e;
	return 0;
}

 int ReadSPE(
       std::vector<double> &x, ///< Spherical, system data, x.e is used. 
	  const char * sysfile ///<input: stream                      
	  )
					  {
    std::ifstream DataStream;  
    DataStream.open(sysfile);
    if (!DataStream)   {               // if the file does not exist
      FatalError("ReadSPE: File "<<sysfile<<" is not found \n"); }
      ReadSPE(x, DataStream);
      DataStream.close();
      return 0;
					  }
  
/**
\page p_pandeya Pandeya Transformation

  Fundamentally, the transformation to the particle hole channel utilizes the commutation identities
\f[ a_{1}^{\dagger}a_{2}^{\dagger}a_{4}a_{3}=-a_{1}^{\dagger}a_{4}a_{2}^{\dagger}a_{3}+\delta_{24}a_{1}^{\dagger}a_{3}=a_{2}^{\dagger}a_{4}a_{1}^{\dagger}a_{3}-\delta_{14}a_{2}^{\dagger}a_{3}\,; 
\f]
the permutation \f$ 1\leftrightarrow 2 \f$ shown here amounts to the action of \f$ \hat{\Theta}_{L\, T}^{(12)}\,\f$ . 
For the summation over indices \f$ \Lambda \f$ and \f$ \sigma \f$ in the above Hamiltonian it is convenient to use a recoupling formula and then utilize the orthogonality conditions of the 3-j symbols. The original two body Hamiltonian 
\f[
  H=\sum_{L\Lambda T\tau}\sum_{{\bf 1}{\bf 2}{\bf 3}{\bf 4}}\, \frac{(1+\delta_{{\bf 1}{\bf 2}})(1+\delta_{{\bf 3}{\bf 4}})}{4}\, V_{L\, T}^{({\bf 1}{\bf 2};{\bf 3}{\bf 4})}\,\left(P_{L\,\Lambda,\, T\,\tau}^{({\bf 1}\,{\bf 2})}\right)^{\dagger}\, P_{L\,\Lambda,\, T\,\tau}^{({\bf 3}\,{\bf 4})}
\f] 
can then be expressed as 
\f[ 
  H=-\frac{1}{4}\sum_{{\bf 1}{\bf 2}{\bf 3}}\,{\delta_{j_{1}j_{3}}\,\sqrt{(1+\delta_{{\bf 1}{\bf 2}})(1+\delta_{{\bf 3}{\bf 2}})}\,\sqrt{2(2j_{2}+1)}}\,\tilde{V}_{0\,0}^{({\bf 2}{\bf 2};{\bf 3}{\bf 1})}\left({\cal M}_{0\,0,\,0\,0}^{({\bf 1}\,{\bf 3})}\right)^{\dagger}\, -
  \frac{1}{4}\sum_{{\bf 1}{\bf 2}{\bf 3}{\bf 4}}\sum_{K\kappa S\sigma}{(2K+1)(2S+1)\sqrt{(1+\delta_{{\bf 1}{\bf 2}})(1+\delta_{{\bf 3}{\bf 4}})}{\tilde{V}}_{KS}^{({\bf 1}{\bf 4};{\bf 3}{\bf 2})}}\,\left({\cal M}_{K\,\kappa,\, S\,\sigma}^{({\bf 1}\,{\bf 4})}\right)^{\dagger}{\cal M}_{K\,\kappa,\, S\,\sigma}^{({\bf 3}\,{\bf 2})}\,, 
\f]
where 
\f[ 
{\tilde{V}}_{KS}^{({\bf 1}{\bf 4};{\bf 3}{\bf 2})}=\sum_{LT}(2L+1)(2T+1) \,\left\{ \begin{array}{ccc}
  1/2 & 1/2 & T\\
  1/2 & 1/2 & S\end{array}\right\} \left\{ \begin{array}{ccc}
  j_{1} & j_{2} & L\\
  j_{3} & j_{4} & K\end{array}\right\} \, V_{L\, T}^{({\bf 1}{\bf 2};{\bf 3}{\bf 4})}\,, 
\f]
transforms the interaction \f$ V\rightarrow\tilde{V}, \f$  from particle-particle to particle-hole channel. The case \f$ K=S=0 \f$ leads to a useful formula for the monopole part 
\f[
{\tilde{V}}_{00}^{({\bf 2}{\bf 4};{\bf 3}{\bf 1})}=-\delta_{j_{1}j_{3}}\,\delta_{j_{2}j_{4}}\,\sum_{LT}\frac{(2L+1)(2T+1)}{2\sqrt{(2j_{1}+1)(2j_{2}+1)}}\, V_{L\, T}^{({\bf 1}{\bf 2};{\bf 3}{\bf 4})}\,. 
\f]
The inverse transformation is obtained using the orthogonality properties of the 6j-symbols 
\f[
  V_{L\, T}^{({\bf 1}{\bf 2};{\bf 3}{\bf 4})}=\sum_{KS}\,(2K+1)(2S+1)\left\{ \begin{array}{ccc}
  1/2 & 1/2 & S\\
  1/2 & 1/2 & T\end{array}\right\} \left\{ \begin{array}{ccc}
  j_{3} & j_{2} & K\\
  j_{1} & j_{4} & L\end{array}\right\} {\tilde{V}}_{KS}^{({\bf 1}{\bf 4};{\bf 3}{\bf 2})}\,. 
\f]
The number of the independent physical constants defining the two-body interactions must be the same in the p-p and p-h channels, meaning that physical limitations similar to, that are directly related to fermionic symmetry, exist for the set of \f$ \tilde{V}\,.\f$ The reverse transformation should lead only to physical \f$ V \f$. 
This is equivalent to a constraint on \f$ \tilde{V} \f$ 
\f[
  {\tilde{V}}_{KS}^{({\bf 1}{\bf 4};{\bf 3}{\bf 2})}=\sum_{K'S'}(-)^{j_{1}+j_{2}-K-K'-S-S'} \,(2K'+1)(2S'+1)\,\left\{ \begin{array}{ccc}
  j_{1} & j_{3} & K'\\
  j_{2} & j_{4} & K\end{array}\right\} \left\{ \begin{array}{ccc}
  1/2 & 1/2 & S'\\
  1/2 & 1/2 & S\end{array}\right\} {\tilde{V}}_{K'S'}^{({\bf 2}{\bf 4};{\bf 3}{\bf 1})}\,. 
\f]
The case of \f$ K=0,\,\, S=0 \f$ is again of particular interest, 
\f[
  {\tilde{V}}_{00}^{({\bf 1}{\bf 4};{\bf 3}{\bf 2})}=-\delta_{j_{1}j_{4}}\,\delta_{j_{2}j_{3}}\,\sum_{KS}\,\frac{(2K+1)(2S+1){\tilde{V}}_{KS}^{({\bf 2}{\bf 4};{\bf 3}{\bf 1})}}{2\sqrt{(2j_{1}+1)(2j_{3}+1)}}\,.
\f]
*/
}
#endif // __INT_CXX__ //

