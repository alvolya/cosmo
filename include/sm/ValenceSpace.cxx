/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file ValenceSpace.cxx
\brief Definitions of valence space and model space for spherical shell model
*/
/*
10/27/09 default read N is changed to particles, parity name change bug fix
*/
#ifndef __VALENCESPACE__CXX__
#define __VALENCESPACE__CXX__
#include <vector>
using std::vector;
#include "QN.cxx"
#include <iostream>
#include <av/iofile.cxx>
#include <av/TextAnalysis.cxx>
#include <sm/label.cxx>
/*! \brief Read data array labeled with "tag" from text file.
*/
template <class v_data>
bool ReadTextData(
v_data *mydata, ///<pointer to a generic output string 
int n, ///<dimension of data 
const char *tag, ///< tag indicating data start
std::ifstream &DataStream ///<stream
){
char junkchar[30];
DataStream.seekg(std::ios::beg); //this line is important, 
    long FPosition=av::FileStringPositionFind(tag, DataStream);
    if (FPosition==-1) return true; //true if fail
    else {
    DataStream>>junkchar;
	for (int i=0;i<n;i++) DataStream>>mydata[i];
	 if (!DataStream.good()) {
        if (DataStream.fail()) FatalError("ReadTextData: Loss of stream integrity while reading \""<<tag<<"\" ");
        if (DataStream.eof()) FatalError("ReadTextData: End-of-file reached while reading \""<<tag<<"\" ");
      }
    return false;
    }
}
/*! \brief Read single data labeled with "tag" from text file.
*/
template <class v_data>
bool ReadTextData(
v_data &mydata, ///<pointer to a generic output string 
const char *tag, ///< tag indicating data start
std::ifstream &DataStream ///<stream
){
return ReadTextData(&mydata,1,tag,DataStream);
}






namespace SM{
/*! \typedef ValenceSpace
\brief Define ValenceSpace: Array of s.p. levels with quantum numbers
This is an array of SPS, each element corresponds to an m-scheme 
position, and contains all quantum numbers related to this position
*/
typedef vector<NonAbelianQN> ValenceSpace;

/*! \typedef ModelSpace
\brief Define ModelSpace: ValenceSpace+symmetry
This is a full ModelSpce, previosly known as Space
*/
typedef SymmetryOperator<ValenceSpace> ModelSpace;

//-----------------Work On Valence Space -----------------------------------------------------
int ReadValenceSpace(ValenceSpace &x, std::ifstream &DataStream){
   int levels;
   if (ReadTextData(levels,"type",DataStream)) FatalError("We need to know the number of s.p. levels");
    x.resize(levels);
    int *j=new int [levels]; 
	//cerr<<levels<<endl;
	if (ReadTextData(j,levels,"spin",DataStream)) {FatalError("We need to know spin j");}
	else for (int i=0;i<levels;i++) x[i].J=j[i]-1;
	if (ReadTextData(j,levels,"isos",DataStream)) {FatalError("We need to know isospin t");}
	else for (int i=0;i<levels;i++) x[i].T=j[i]-1;
	if (ReadTextData(j,levels,"mqn",DataStream)) 
	{cerr<<"Main q.n. are unknown\n"; for (int i=0;i<levels;i++) x[i].N=0;} 
	else {for (int i=0;i<levels;i++) x[i].N=j[i];}
	
	if (ReadTextData(j,levels,"orb",DataStream)) {cerr<<"orbital q.n. are unknown\n";
	for (int i=0;i<levels;i++) x[i].L=0;
	}
	else for (int i=0;i<levels;i++) x[i].L=j[i];
	for (int i=0;i<levels;i++) x[i].tag=SPLabel(x[i].N, x[i].L, x[i].J);
	delete [] j;
	return 0;
}

int ReadValenceSpace(ValenceSpace &x, const char* myfile){
  std::ifstream DataStream;
  DataStream.open(av::FileName(myfile,"SMINTPATH",".smin").c_str(), std::ios::in|std::ios::binary);
   INF(std::cerr<<"Reading Valence space from: "<<av::FileName(myfile,"SMINTPATH",".smin")<<"\n");
   //DataStream.open(myfile, std::ios::in|std::ios::binary);
   ReadValenceSpace(x,DataStream);
   DataStream.close();
   return 0;
}

int WriteValenceSpace(ValenceSpace &x, std::ofstream &DataStream){
    DataStream<<"! Spins and Isospins (with capital) correspond to actual values (not 2j+1 as in the input) \n";
    DataStream<<"types "<<x.size()<<"\n";
	DataStream<<"Spins ";
	for (int i=0;i<x.size();i++) DataStream<<x[i].J<<" "; DataStream<<"\n";
	DataStream<<"Isosp "; 
	for (int i=0;i<x.size();i++) DataStream<<x[i].T<<" "; DataStream<<"\n";
	DataStream<<"mqn "; 
	for (int i=0;i<x.size();i++) DataStream<<x[i].N<<" "; DataStream<<"\n";
	DataStream<<"orbital "; 
	for (int i=0;i<x.size();i++) DataStream<<x[i].L<<" "; DataStream<<"\n";
	return 0;
}





int WriteValenceSpace(ValenceSpace &x, const char* myfile){
  std::ofstream DataStream;
  DataStream.open(myfile);
   WriteValenceSpace(x,DataStream);
   DataStream.close();
   return 0;
}

//Work on Model Space ----------------------------------------------------------------
std::string LabelModelSpace(ModelSpace &x) {
if (x.size()==0) FatalError("The valence space is empty");
if (x[0].J&1) x.tag="f"; else x.tag="b";
char  stmp[256];
	for (int i=0;i<x.size();i++) {
	std::snprintf (stmp,256, "l%i", x[i].J+1);
	 x.tag+=stmp;
  }
   
	std::snprintf (stmp,256, "M%iT%iN%i", int(x.Jz),int(x.Tz),x.N);
	x.tag+=stmp;
	if (x.P) x.tag+="-"; else x.tag+="+";
    return x.tag;
}

int ReadModelSpace(ModelSpace &x, std::ifstream &DataStream){
	ReadValenceSpace(x,DataStream);
   if (ReadTextData(x.N,"particles",DataStream)) {
   cerr<<"Particle number is not defined, enter N= ";
    std::cin>>x.N;
   }
   
   if (ReadTextData(x.Jz,"JZ",DataStream)) {
   cerr<<"Spin projection is not defined, enter Jz= ";
       std::cin>>x.Jz;
   }
   
   if (ReadTextData(x.Tz,"TZ",DataStream)) {
   cerr<<"Isopin projection is not defined, enter Tz= ";
       std::cin>>x.Tz;
   }
   
   if (ReadTextData(x.P,"parity",DataStream)) {
   cerr<<"Parity is not defined, enter (+1 or -1) P= ";
    std::cin>>x.P;
    }
	if (x.P>=0) x.P=0; else x.P=1; //read parity in human format
	
	 if (ReadTextData(x.tag,"label",DataStream)) {
   cerr<<"label is not defined; new label= ";
   cerr<<LabelModelSpace(x)<<endl;
    }
	

//	for (int i=0;i<levels;i++) x[i].tag=SPLabel(x[i].N, x[i].L, x[i].J);

	return 0;
}
int ReadModelSpace(ModelSpace &x, const char* myfile){
  std::ifstream DataStream;
  DataStream.open(av::FileName(myfile,"SMINTPATH",".smin").c_str(), std::ios::in|std::ios::binary);
   INF(std::cerr<<"Reading Model space from: "<<av::FileName(myfile,"SMINTPATH",".smin")<<"\n");
   //DataStream.open(myfile, std::ios::in|std::ios::binary);
   ReadModelSpace(x,DataStream);
   DataStream.close();
   return 0;
}

int WriteModelSpace(ModelSpace &x, std::ofstream &DataStream){
	WriteValenceSpace(x,DataStream);
	DataStream<<"label "<<x.tag<<"\n";
	DataStream<<"particles "<<x.N<<"\n";
	DataStream<<"JZ "<<x.Jz<<"\n";
	DataStream<<"TZ "<<x.Tz<<"\n";
	DataStream<<"parity ";
	if (x.P) DataStream<<"-1 \n"; else DataStream<<"1 \n";
	return 0;
}

int WriteModelSpace(ModelSpace &x, const char* myfile){
  std::ofstream DataStream;
  DataStream.open(myfile);
   WriteModelSpace(x,DataStream);
   DataStream.close();
   return 0;
}

/*! \brief makes an oscillator valence space with all states up to Nmax quanta.
assuming all particles have spin s=s2/2 and isospin t2/2
 */
int MakeNmaxValenceSpace(
    SM::ValenceSpace& x, ///< output valence space
    int Nmax, ///< input Nmax integer
    int s2, ///< input spin of each particle
    int t2 ///< input isospin of each particle
    )
{
    SM::NonAbelianQN qtmp;
    int l(0);
    //double j(0.0);
    for (int i = 0; i <= Nmax; i++)//i is the number of quanta in the shell.
    {
        for (int n = i/2; n >= 0; n--)//number of nodes.
        {
            l = i - 2 * n;
           
            qtmp.N = n;
            qtmp.L = l;
            qtmp.T = t2/2.0; //can use spin variable
            qtmp.S = s2/2.0;
            
            for (int j2=std::abs(2*l-s2); j2<=(2*l+s2); j2+=2) {
             qtmp.J = j2/2.0;
            x.push_back(qtmp);
            }
        }
    }
return x.size();
}

int MakeNmaxValenceSpace(SM::ValenceSpace& x,int Nmax)
{
return MakeNmaxValenceSpace(x,Nmax,1,1); //just return specific case 
}

/*! \brief makes an oscillator valence space with all states up to Nmax quanta.
 */
//int MakeNmaxValenceSpace(SM::ValenceSpace& x,int Nmax)
//{
//    SM::NonAbelianQN qtmp;
//    int l(0);
//    double j(0.0);
//    for (int i = 0; i <= Nmax; i++)//i is the number of quanta in the shell.
//    {
//        for (int n = i/2; n >= 0; n--)//number of nodes.
//        {
//            l = i - 2 * n;
//            j = (2*l-1)/2.;
//            qtmp.N = n;
//            qtmp.L = l;
//            qtmp.J = j;
//            qtmp.T = 0.5;
//            qtmp.S = 0.5;
//            
//            if ( l != 0 ) x.push_back(qtmp);//s-levels don't have j=l-1/2.
//            qtmp.J = (2*l + 1)/2.;
//            x.push_back(qtmp);
//        }
//    }
//return x.size();
//}



//Print valece space
std::ostream& operator << (std::ostream &stream,
                                const SM::ValenceSpace& x)
    {
    for (int i=0;i<x.size();i++)
        stream<<"n="<<x[i].N<<" l="<<x[i].L<<" j="<<x[i].J<<endl;
     return stream;
}

    int FindNmax(const SM::ValenceSpace& x) {
        int Nmax=0; //start with lowest possible 0
        int nn;
        for (int i=0;i<x.size();i++) {
            nn=2*x[i].N+x[i].L;
            if (Nmax<nn) Nmax=nn;
        }
        return Nmax;
        }
    int FindNmin(const SM::ValenceSpace& x) {
        if (x.size()==0) return 0;
        int Nmin=2*x[0].N+x[0].L; //start with first value
        int nn;
        for (int i=1;i<x.size();i++) {
            nn=2*x[i].N+x[i].L;
            if (Nmin>nn) Nmin=nn;
        }
        return Nmin;
    }


}//end SM namespace
#endif




