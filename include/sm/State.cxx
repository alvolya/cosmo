/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


#ifndef __STATE__CXX__
#define __STATE__CXX__
/*!
\author Alexander Volya  <http://www.volya.net>
\file State.cxx
\brief Spherical Shell Model Eigen States database driver
\todo This State requires a lot of revisions 
use memory or use file for a given state (perhaps tensors in files?)
*/
//FEigenstate is not used
/*
03/13/2006  
          change array to a tensor structure
	  reference to mbs is only included if MBS is defined

class State
stream PrintState(stream, State)
stream ReadState(stream, State)
operaoor << State
operator << Tensor<1,State>
EnergyStateRead(energy, tensor, filename, id)
Read   Write  binary
Render Pring  text
04/14/08 introduce EBinary into state, corresponds to energy in binary file

 */
//#include <sm/ValenceSpace.cxx>
#include <cstdlib> 
using std::getenv;
//use getenv to read path of the interaction files
#include <av/iofile.cxx> //reading file with path
#include <av/TextAnalysis.cxx>
#include <sm/QN.cxx>
#include <sm/label.cxx>
#include <tav/tensor.cxx>
#include <rw/basic-rw.cxx> 
//using rw::Read;
using std::strcmp;
namespace SM { 

  //! this is a database structure for eigen states   
 class State: public QN {
    public:   
   char label[256]; //this is the EE file
   int ID; //this is the number of the state in EE file
   double E;
   double EBinary; //! Energy in binary file
#ifdef __MANYBODYSTATES__CXX__
   Many_Body_States *basis; //pointer to basis states
#endif
   char basisfile[256]; //file where basis is stored
   //long int n; //dimension
   //double *z; //the vector itself
   tav::Tensor<1,double> z; //this is vector itself
   // State() {z=NULL;} //default constructor
  // ~State() {if (z!=NULL) delete [] z;} //default destructor
   int Read(char *filename);
   int Write(char *filename);
    };

std::ostream& PrintState(std::ostream& stream, State &EV)
{
    stream.width(10);
    stream<<EV.E<<' ';
    stream.width(6);
    stream<<EV.J<<' ';
    stream.width(6);
    stream<<EV.T<<' ';
    stream.width(6);
    stream<<EV.L<<' ';
    stream.width(6);
    stream<<EV.N<<' ';    
    stream<<EV.label<<'\t';
    stream<<EV.ID<<' ';     
    //  stream<<EV.N<<'\t';
  //  stream<<EV.Jz<<'\t';
  //  stream<<EV.Tz<<'\t';
  //  stream<<EV.P<<'\t';
  //  stream<<EV.E-EV[0].E;
    //  stream<<'\n';
  stream.width();
        return stream;
}


std::ostream& operator << (std::ostream& stream, State &EV)
{
  PrintState(stream, EV);
  stream<<'\n';
  return stream;
}

  //ouptut entire vector technically only needed to show excitation energies
std::ostream& operator << (std::ostream& stream, tav::Tensor<1,State> &EV)
{
  stream<<"# -0.5 means unknown \n";
  stream<<"#ENERGY  J  T L N datafile location excitation "<<endl;
  stream<<EV.dim<<'\n';
  for (int i=0;i<EV.dim;i++) {
    PrintState(stream, EV[i]);
    stream<<EV[i].E-EV[0].E;
  stream<<'\n';
   }
        return stream;
}

 std::istream& RenderState (std::istream& stream, State &EV)
{
     stream >> EV.E;
    double rspin;
      stream>>EV.J;
      stream>>EV.T;
      stream>>EV.L;
      stream>>EV.N;
      stream >> EV.label;
      ReadQNLabel(EV.label, EV);
      stream >> EV.ID;
      //read QN
      // inf>>EV.N;
      //  inf>>EV.Tz;
      //  inf>>EV.P;
      //end reading
      //inf >> EV;
      //stream.getline(temporary_name,256);
        return stream;
}

  std::istream& operator >> (std::istream& stream, State &EV)
{
  RenderState(stream, EV);
  stream.getline(temporary_name,256);
        return stream;
}



int EnergyStateRead(double &E, tav::Tensor<1,double> &EV, const char *filename, int ID) {
  std::ifstream inf;
      inf.open(filename ,std::ios::binary); //open binary
      int currentID;
    rw::Read(currentID,inf); //read number of levels
      if ((currentID-1)<ID) 
	cerr<<"We have a problem: state is not in file EE"<<endl;
      for (int j=0;j<=ID;j++) {//sequentially read file
      double ErrorE;
      double ErrorV;
    rw::Read(E,inf);
    rw::Read(ErrorE,inf);
    rw::Read(ErrorV,inf);
	EV.Read(inf);	
      
      }
      inf.close();
      return ID;
}

  int Read(tav::Tensor<1,State> &EV, const char *filename) 
  {
  std::ifstream inf;
  //inf.open(SM::AddExtension(filename,".ev"));
  inf.open(av::FileName(filename,"SMPATH",".ev").c_str());
 // INF(std::cerr<<"Reading Interactions from \n"<<av::FileName(sysfile,"SMINTPATH",".INT")<<std::endl);
  int n; //number of states stored in file
  //inf || n;
  av::SkipComments(inf);
  inf>>n;
  //cerr<<n<<endl;
  //Pause();
  EV.Resize(n);
  for (int i=0;i<n;i++) {
    RenderState(inf, EV[i]);
    inf.getline(temporary_name,256); //ignore the rest of the line
  } 
  inf.close();
  //we now need to read eigenvectors sequential read
  char currentname[256];
  int currentID=-1;
  strcpy(currentname, EV[0].label); //first file to open
  for (int i=0;i<n;i++) {
    if ((inf.is_open())&&(strcmp(currentname,EV[i].label))&&(currentID==EV[i].ID+1))
      {
          rw::Read(EV[i].EBinary,inf);
	EV[i].z.Read(inf);
      } //this is if all is good and running
    else {
      if ((inf.is_open())) inf.close(); //close file first
      strcpy(currentname, EV[i].label); //copy label
      inf.open(SM::AddExtension(currentname,".HH.EE"),std::ios::binary); //open binary
        rw::Read(currentID,inf); //read number of levels
      if ((currentID-1)<EV[i].ID) 
      {FatalError("Read(*.ev): We have a problem: state is not in file EE"); }
      for (int j=0;j<=EV[i].ID;j++) {//sequentially read file
          rw::Read(EV[i].EBinary,inf);
	EV[i].z.Read(inf);	
      }
    }//end if 
    currentID=EV[i].ID; //fix ID, label should be correct at this point  
    if (fabs(EV[i].EBinary-EV[i].E)>1E-3) 
      WMessage("The energy of state "<<i<<" in file "<<EV[i].label<<".ev \n is changed, will keep updated value.");
  }
  return 1;
  }
  /*! \Read ev file and binary eigenvalues, create array of states with energies as corrections
  */
  int ReadCorrections(tav::Tensor<1,State> &EV, char *filename) 
  {
    std::ifstream inf;
    inf.open(SM::AddExtension(filename,".ev"));
    int n; //number of states stored in file
    //n=av::LineCount(SM::AddExtension(filename,".ech"));
    av::SkipComments(inf);
    inf>>n;
    //inf || n;
    EV.Resize(n);
    for (int i=0;i<n;i++) {
      RenderState(inf, EV[i]);
      inf.getline(temporary_name,256); //ignore the rest of the line
    } 
    inf.close();
  //we now need to read eigenvectors sequential read
    char currentname[256];
    int currentID=-1;
    double currentE;
    strcpy(currentname, EV[0].label); //first file to open
    for (int i=0;i<n;i++) {
      if ((inf.is_open())&&(strcmp(currentname,EV[i].label))&&(currentID==EV[i].ID+1))
      {
          rw::Read(currentE,inf);
        EV[i].z.Read(inf);
      } //this is if all is good and running
      else {
        if ((inf.is_open())) inf.close(); //close file first
        strcpy(currentname, EV[i].label); //copy label
        inf.open(SM::AddExtension(currentname,".HH.EE"),std::ios::binary); //open binary
          rw::Read(currentID,inf); //read number of levels
        if ((currentID-1)<EV[i].ID) 
        {FatalError("ReadCorrection(*.ev): We have a problem: state "<<currentID<<" is not in file EE"); }
        for (int j=0;j<=EV[i].ID;j++) {//sequentially read file
            rw::Read(currentE,inf);
          EV[i].z.Read(inf);	
        }
      }//end if 
      currentID=EV[i].ID; //fix ID, label should be correct at this point  
      if (fabs(currentE-EV[i].E)>1E-3) {
        //CriticalStop("The energy of state "<<i<<" in file "<<EV[i].label<<".ev \n is changed, will keep updated value.");
        //currentE-entry is bindary entry file, EV[i].E entry in database
        EV[i].E-=currentE;
      //cerr<<"Use updated energy value"<<endl;
      }
      else EV[i].E=0.0;
    }

    return 1;
  }
  
  /*
  template <class cnumber>
  int Read(cnumber &T, std::ifstream &inf) {
  inf.read((char*)(&T),sizeof(cnumber)); 
  return 1;
  }
  */
  
}

#endif //__STATE__CXX__

