/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>

\file ReadSPInt.cxx
\brief Read spherical interactions and convert them to single particle format on-fly
\ingroup gp_nr
*/
#ifndef __READSPINT__CXX__
#define __READSPINT__CXX__
#include <sm/spherical2sp.cxx>
  /*! \brief Read interaction in spherical format and output single-particle format, 
  input file includes scalings in vector of scalings and isospins
   */
  int ReadSPINT(
tav::STensor<spsint,double> &VSP,  ///< Output map of elements V[1234]  (12) (34)-ordered pairs
SM::SingleParticleStates & sps, ///< s.p. states matrix
const std::vector<double> &scalings, //scalings for columbs in file
const std::vector<int> &ispins, //isospins 1-NN, -1 PP, 0 PN, else isoinvariant
std::ifstream &DataStream ///<input: stream
                      ) 
  {
    const int maxid=CountLevels(sps);
    const double minread=1E-8; //minimum value of matrix element to read
    const int MaxLineLength=1024;
    char   strFileLine[MaxLineLength]; // One line buffer...
    long FPosition;
//    DataStream.seekg(std::ios::beg); //this line is important
    FPosition=av::SkipComments(DataStream);
    if (FPosition==-1) FatalError("No data is available in stream");    
    int intlines=0;
    DataStream>>intlines;  //in testing version skip the rest
    DataStream.getline(strFileLine, MaxLineLength, '\n');
    intlines=abs(intlines); //just in case this is oxbash file
    int ijk[6]; 
    double htmp;
    double melement;
    bool rejectline=false;
    for (int il=0;il<intlines;il++) {
    rejectline=false;
      DataStream>>ijk[2];
      DataStream>>ijk[3];
      DataStream>>ijk[0];
      DataStream>>ijk[1];
      DataStream>>ijk[4];
      DataStream>>ijk[5];
  //DataStream>>VX;
  //convert to from oxbash notation
      if ((--ijk[0]>=maxid)||(--ijk[1]>=maxid)||(--ijk[2]>=maxid)||(--ijk[3]>=maxid)) rejectline=true;
//      --ijk[1];
//      --ijk[2];
//      --ijk[3];
  // Note OXBASH ORDERING IS FROM HIGHER TO LOWER
 // if (ijk[0]>ijk[1]) tav::Swap(ijk[0],ijk[1]);
 // if (ijk[2]>ijk[3]) tav::Swap(ijk[2],ijk[3]);
 // if (ijk[0]>ijk[2]) {tav::Swap(ijk[0],ijk[2]); tav::Swap(ijk[1],ijk[1]);
 // we are now ready to read data
    for (int ic=0;ic<scalings.size();ic++)
    {
    DataStream>>htmp;
    htmp*=scalings[ic];
    if (rejectline) continue; //bad line after read continue
    if (std::abs(htmp)<minread) continue; //matrix element is too small
    if (std::abs(ispins[ic])<2) AddPPVsp(VSP, sps, ijk, htmp,ispins[ic]);
    else AddPPVsp(VSP, sps, ijk, htmp);
    }
    DataStream.getline(strFileLine, MaxLineLength, '\n'); //this is used to ignore the rest of data in line
    CheckInputStream(DataStream);
    }
    return intlines;
  };
/** \brief Read Oxbash-format interaction from multi-column file
 */

int ReadSPINT(
tav::STensor<spsint,double> &VSP,  ///< Output map of elements V[1234]  (12) (34)-ordered pairs
SM::SingleParticleStates & sps, ///< s.p. states matrix
const std::vector<double> &scalings, //scalings for columbs in file
const std::vector<int> &ispins, //isospins 1-NN, -1 PP, 0 PN, else isoinvariant
  const char *sysfile ///<file name
                      )
  {

    std::ifstream DataStream;  
    DataStream.open(sysfile);
    if (!DataStream)   {               // if the file does not exist
      FatalError("ReadSPINT: File "<<sysfile<<" is not found \n"); }
      ReadSPINT(VSP, sps, scalings, ispins, DataStream);
      DataStream.close();
      return 0;
  }
#endif //__READSPINT__CXX__
