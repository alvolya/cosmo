/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file MakeMBS.cxx
\brief In this file we use various algorithms to make the Many-Body States
\ingroup gp_mbs
*/
#ifndef __MAKEMBS__CXX__
#define __MAKEMBS__CXX__
//#include <tav/sort.cxx>
//using tav::QuickSort;
//#include <tav/sort.cxx>
//using tav::QuickSort;
#include <tav/pconfiguration.cxx>
using tav::FractionalConfigurationPosition;
namespace SM {

/*! \brief Make many-body states, bosons-fermions
Works for bosons and fermions, abelian addition implemented as particles are moved
\ingroup gp_mbs
*/
uint_nbasis MakeMBS(const char *filename, int N, int Jz, int Tz, int P, SingleParticleStates &sps)
{
  uint_nbasis n=0;
  int type=sps[0].J%2;
 ofstream outf;
  outf.open(filename, ios::binary);
  outf.write((char*)(&N),sizeof(int));
  outf.write((char*)(&n),sizeof(uint_nbasis)); //write zero for now

  spin TZ;
  spin JZ;
  bool LL; //abelian spin and parity variables

  
  spsint *currentstate=new spsint [N]; 
  //this is my current state 

  //fermions
  if (type==1) {
    //make an initial distribution
    if (sps.size()<N) {cerr<<"fermions do not fit in this space"<<endl; return 0;}
    n=0;
    JZ=0;TZ=0;LL=true;
    for (int i=0;i<N;i++) {
      currentstate[i]=i;
      JZ+=sps[i].Jz;
      TZ+=sps[i].Tz;
      LL^=sps[i].P;}
    
    // cout<<"initial state is done"<<endl;
    //[xxxx-------] initial distribution
    
    //do an increment
    bool OKflag=true;
    
    while (OKflag) {
      //control conservation laws
      if ((JZ==Jz)&&(TZ==Tz)&&(LL^P)) 
     {
       //here the wanted state is ready
       //Numbering starts from 0
       //  cout<<n<<" JZ="<<JZ<<" TZ="<<TZ<<" L="<<LL<<" : ";
       //for (int i=0;i<N;i++) cout<<int(currentstate[i])<<" ";
       //cout<<endl;
       outf.write((char*)(currentstate), N*sizeof(spsint));
	 n++;
	 // if ((n%1000)==0) cout<<n<<endl;
    //n gives number of states
}
   

      OKflag=false;
    for (int i=N-1;i>=0;i--) { 
      //go in reversed order
      //  if (currentstate[i]+1==sps.size()) continue; //we are at end [x-xxx----X]
      // if ((i!=(N-1))&&((currentstate[i]+1)==currentstate[i+1])) continue;
      // we need to move i-th particle and dump all others after it
      // [x--xX---xxx]->[x--x-Xxxx---]
      //number of spaces I have after it is sps.size()-1-currentstate[i]
      //number of particle including X is N-i
      if ((sps.size()-1-currentstate[i])<(N-i)) continue; //we do not fit;
      //there is a trick here we need to change current i last
      for (int j=N-1;j>=i;j--) {
	//remove all quantum numbers
	JZ-=sps[currentstate[j]].Jz;
	TZ-=sps[currentstate[j]].Tz;
	LL^=sps[currentstate[j]].P;
	currentstate[j]=currentstate[i]+j+1-i; //move particle
	//put quantum numbers back
	JZ+=sps[currentstate[j]].Jz;
	TZ+=sps[currentstate[j]].Tz;
	LL^=sps[currentstate[j]].P;

}
      // n++;


      OKflag=true;
      break;
    }
    }; 

  }
 //bosons
  if (type==0) {
    //make an initial distribution
    n=0;
  JZ=0;TZ=0;LL=true;
    for (int i=0;i<N;i++) {
      currentstate[i]=0;
      JZ+=sps[0].Jz;
      TZ+=sps[0].Tz;
      LL^=sps[0].P;
    }  
    // cout<<"initial state is done"<<endl;
    //[xxxx-------] initial distribution
    
    //do an increment
    bool OKflag=true;
    
    while (OKflag) {

     //control conservation laws
      if ((JZ==Jz)&&(TZ==Tz)&&(LL^P)) 
     {
       //here the wanted state is ready
       //Numbering starts from 0
       //  cout<<n<<" JZ="<<JZ<<" TZ="<<TZ<<" L="<<LL<<" : ";
       //for (int i=0;i<N;i++) cout<<int(currentstate[i])<<" ";
       //cout<<endl;
       outf.write((char*)(currentstate), N*sizeof(spsint));
    n++;
    //n gives number of states
}
   
 
      OKflag=false;
    for (int i=N-1;i>=0;i--) { //go in reversed order
      //  if (currentstate[i]+1==sps.size()) continue; //we are at end [x-xxx----X]
      // if ((i!=(N-1))&&((currentstate[i]+1)==currentstate[i+1])) continue;
      // we need to move i-th particle and dump all others after it
      // [x--xX---xxx]->[x--x-Xxxx---]
      //number of spaces I have after it is sps.size()-1-currentstate[i]
      //number of particle including X is N-i
      if (currentstate[i]==(sps.size()-1)) continue; //we are at end;
      //there is a trick here we need to change current i last
      for (int j=N-1;j>=i;j--)  {
	JZ-=sps[currentstate[j]].Jz;
	TZ-=sps[currentstate[j]].Tz;
	LL^=sps[currentstate[j]].P;
	 //move particle
	currentstate[j]=currentstate[i]+1;
	//put quantum numbers back
	JZ+=sps[currentstate[j]].Jz;
	TZ+=sps[currentstate[j]].Tz;
	LL^=sps[currentstate[j]].P;


      } 
      OKflag=true;
      break;
    }
    }; 

  } //---------------END BOSONS-------------------
  delete [] currentstate;
  //now I need to overwrite n in file
  outf.seekp(0); //go to the beguinning of the file
 outf.write((char*)(&N),sizeof(int));
  outf.write((char*)(&n),sizeof(uint_nbasis)); //overide new n and N
outf.close();
  return n;
}

/*! \brief Make Many-Body states and save into file, fast method.
This is a fast method because particle distribution search is terminated subset has no possibility to satisfy abelian QN
\ingroup gp_mbs
*/
uint_nbasis MakeMBS2(
const char *filename, ///<output file
int N, int Jz, int Tz, int P, ///<system 
SingleParticleStates &sps ///<single-particle states
)
{
    //process status variable if defined
#ifdef CSTATUS
    av::status S_MBS("MBS");
#ifdef FSTATUS
    S_MBS.SetFile();
#endif
#endif
    //-----------end status
  uint_nbasis n=0;
  int type=sps[0].J%2;
 ofstream outf;
  outf.open(filename, ios::binary);
  outf.write((char*)(&N),sizeof(int));
  outf.write((char*)(&n),sizeof(uint_nbasis)); //write zero for now

  spin TZ;
  spin JZ;
  bool LL; //abelian spin and parity variables

  //create list Mmax
  //I need a table given that ni particles what maximum M and T they could have
  int *Mmax= new int [N+1];
  Mmax[0]=0; //no particles no spin
  //I have an array sps[i].Jz at my disposal
  int *mlist=new int [sps.size()];
  for (int i=0;i<sps.size();i++) mlist[i]=sps[i].Jz; //copy list of M
  //QuickSort(sps.size(),mlist); //we sorted the list, print it
  std::sort(mlist, mlist + sps.size());
  //on one level mlist will be -j -j+1 ... j 
  //   for (int i=0;i<sps.size();i++) cout<<i<<" "<<mlist[i]<<endl;
  for (int in=1;in<=N;in++) Mmax[in]=Mmax[in-1]+mlist[sps.size()-in];
  // delete [] mlist; //finish up and delete sorted list of spins
  //  for (int in=0;in<=N;in++) cout<<in<<" ["<<Mmin[in]<<", "<<Mmax[in]<<"] "<<endl;

 int *Tmax= new int [N+1];
  Tmax[0]=0; //no particles no spin
  //I have an array sps[i].Jz at my disposal
  // int *mlist=new int [sps.size()];
  for (int i=0;i<sps.size();i++) mlist[i]=sps[i].Tz; //copy list of M
  //QuickSort(sps.size(),mlist); //we sorted the list, print it
  std::sort(mlist, mlist + sps.size());
  //on one level mlist will be -j -j+1 ... j 
  //   for (int i=0;i<sps.size();i++) cout<<i<<" "<<mlist[i]<<endl;
  for (int in=1;in<=N;in++) Tmax[in]=Tmax[in-1]+mlist[sps.size()-in];
  delete [] mlist; //finish up and delete sorted list of spins
  //  for (int in=0;in<=N;in++) cout<<in<<" ["<<Mmin[in]<<", "<<Mmax[in]<<"] "<<endl;


  spsint *currentstate=new spsint [N]; 
 
  //fermions
  if (type==1) {
    //make an initial distribution
    if (sps.size()<N) {cerr<<"fermions do not fit in this space"<<endl; return 0;}
 //this is my current state 
  int *MSUM=new int [N];
  int *TSUM=new int [N];
    n=0;
    JZ=0;TZ=0;LL=true;
    for (int i=0;i<N;i++) {
      currentstate[i]=i;
      JZ+=sps[i].Jz;
      TZ+=sps[i].Tz;
      LL^=sps[i].P;}
    
    // cout<<"initial state is done"<<endl;
    //[xxxx-------] initial distribution
    
    //do an increment
    bool OKflag=true;
    
    while (OKflag) {
      //control conservation laws
      if ((JZ==Jz)&&(TZ==Tz)&&(LL^P)) 
     {
       //here the wanted state is ready
       //Numbering starts from 0
       //  cout<<n<<" JZ="<<JZ<<" TZ="<<TZ<<" L="<<LL<<" : ";
       //for (int i=0;i<N;i++) cout<<int(currentstate[i])<<" ";
       //cout<<endl;
       outf.write((char*)(currentstate), N*sizeof(spsint));
	 n++;
     //cerr<<n<<": "<<mbs_id(N, currentstate, sps.size())<<" : "<<int(currentstate[0])<<" "<<int(currentstate[1])<<" "<<sps.size()<<" "<<fmbs_id(N, currentstate, sps.size())<<endl;
#ifdef CSTATUS
         if (S_MBS()) {S_MBS.stream()<<" n=";
         av::PrintInt(S_MBS.stream(),n);
         S_MBS.stream()<<" "<<tav::FractionalConfigurationPosition(sps.size(),N,currentstate)*100.<<"%";
         }
        #endif
	// if ((n%1000000)==0) cout<<" "<<n<<" "<<int(currentstate[0])<<"/"<<sps.size()<<endl;
    //n gives number of states
}
   

      OKflag=false;
      //prepare MSUM;
      MSUM[0]=0;
    for (int in=0;in<N-1;in++) MSUM[in+1]=MSUM[in]+sps[currentstate[in]].Jz;
    TSUM[0]=0;
    for (int in=0;in<N-1;in++) TSUM[in+1]=TSUM[in]+sps[currentstate[in]].Tz;
    //for (int in=0;in<N;in++) cout<<"-------------- "<<MSUM[in]<<endl;
    for (int i=N-1;i>=0;i--) {
  //go in reversed order
      /* I have N-i particles moving starting from I 
      [x-xx-x-Xxx-x-x]
       0 1    i     N-1                 
      say MSUM[i]=sps[0].Jz+sps[1].Jz+..sps[i-1].Jz
      it only worth to proceed if Mmin[N-i]<=Jz-MSUM[i]<=Mmax[N-i]
      */
      
      if (abs(Jz-MSUM[i])>Mmax[N-i]) {continue;}
      if (abs(Tz-TSUM[i])>Tmax[N-i]) {continue;}
      // if (abs(Jz+MSUM[i])>Mmax[N-i]) continue;
      
    
      //  if (currentstate[i]+1==sps.size()) continue; //we are at end [x-xxx----X]
      // if ((i!=(N-1))&&((currentstate[i]+1)==currentstate[i+1])) continue;
      // we need to move i-th particle and dump all others after it
      // [x--xX---xxx]->[x--x-Xxxx---]
      //number of spaces I have after it is sps.size()-1-currentstate[i]
      //number of particle including X is N-i
      if ((sps.size()-1-currentstate[i])<(N-i)) continue; //we do not fit;
      //there is a trick here we need to change current i last
      for (int j=N-1;j>=i;j--) {
	//
	//remove all quantum numbers
	JZ-=sps[currentstate[j]].Jz;
	TZ-=sps[currentstate[j]].Tz;
	LL^=sps[currentstate[j]].P;
	currentstate[j]=currentstate[i]+j+1-i; //move particle
	//put quantum numbers back
	JZ+=sps[currentstate[j]].Jz;
	TZ+=sps[currentstate[j]].Tz;
	LL^=sps[currentstate[j]].P;

}
      // n++;


      OKflag=true;
      break;
    }
    }; 
    delete [] MSUM;
    delete [] TSUM;
  }
 //bosons
  if (type==0) {
    //make an initial distribution
    n=0;
  JZ=0;TZ=0;LL=true;
    for (int i=0;i<N;i++) {
      currentstate[i]=0;
      JZ+=sps[0].Jz;
      TZ+=sps[0].Tz;
      LL^=sps[0].P;
    }  
    // cout<<"initial state is done"<<endl;
    //[xxxx-------] initial distribution
    
    //do an increment
    bool OKflag=true;
    
    while (OKflag) {

     //control conservation laws
      if ((JZ==Jz)&&(TZ==Tz)&&(LL^P)) 
     {
       //here the wanted state is ready
       //Numbering starts from 0
       //  cout<<n<<" JZ="<<JZ<<" TZ="<<TZ<<" L="<<LL<<" : ";
       //for (int i=0;i<N;i++) cout<<int(currentstate[i])<<" ";
       //cout<<endl;
       outf.write((char*)(currentstate), N*sizeof(spsint));
    n++;
    //n gives number of states
}
   
 
      OKflag=false;
    for (int i=N-1;i>=0;i--) { //go in reversed order
      //  if (currentstate[i]+1==sps.size()) continue; //we are at end [x-xxx----X]
      // if ((i!=(N-1))&&((currentstate[i]+1)==currentstate[i+1])) continue;
      // we need to move i-th particle and dump all others after it
      // [x--xX---xxx]->[x--x-Xxxx---]
      //number of spaces I have after it is sps.size()-1-currentstate[i]
      //number of particle including X is N-i
      if (currentstate[i]==(sps.size()-1)) continue; //we are at end;
      //there is a trick here we need to change current i last
      for (int j=N-1;j>=i;j--)  {
	JZ-=sps[currentstate[j]].Jz;
	TZ-=sps[currentstate[j]].Tz;
	LL^=sps[currentstate[j]].P;
	 //move particle
	currentstate[j]=currentstate[i]+1;
	//put quantum numbers back
	JZ+=sps[currentstate[j]].Jz;
	TZ+=sps[currentstate[j]].Tz;
	LL^=sps[currentstate[j]].P;


      } 
      OKflag=true;
      break;
    }
    }; 

  } //---------------END BOSONS-------------------
  delete [] currentstate;
  delete [] Tmax;
  delete [] Mmax;
  //now I need to overwrite n in file
  outf.seekp(0); //go to the beguinning of the file
 outf.write((char*)(&N),sizeof(int));
  outf.write((char*)(&n),sizeof(uint_nbasis)); //overide new n and N
//  cerr<<"writing n="<<n;
outf.close();
  return n;
}

/* We make many-body states with rejection by the number of qunta
*/
uint_nbasis MakeMBS2(
const char *filename, ///<output file
int N, int Jz, int Tz, int P,
int HW, ///< number of oscillator quanta
SingleParticleStates &sps ///<single-particle states
)
{
    //process status variable if defined
#ifdef CSTATUS
    av::status S_MBS("MBS");
#ifdef FSTATUS
    S_MBS.SetFile();
#endif
#endif
    //-----------end status
  uint_nbasis n=0;
  int type=sps[0].J%2;
 ofstream outf;
  outf.open(filename, ios::binary);
  outf.write((char*)(&N),sizeof(int));
  outf.write((char*)(&n),sizeof(uint_nbasis)); //write zero for now

  spin TZ;
  spin JZ;
  bool LL; //abelian spin and parity variables
  int WW; //number of quanta
cerr<<"Number of quanta "<<HW<<endl;
  //create list Mmax
  //I need a table given that ni particles what maximum M and T they could have
  int *Mmax= new int [N+1];
  Mmax[0]=0; //no particles no spin
  //I have an array sps[i].Jz at my disposal
  int *mlist=new int [sps.size()];
  for (int i=0;i<sps.size();i++) mlist[i]=sps[i].Jz; //copy list of M
  //QuickSort(sps.size(),mlist); //we sorted the list, print it
  std::sort(mlist, mlist + sps.size());
  //on one level mlist will be -j -j+1 ... j 
  //   for (int i=0;i<sps.size();i++) cout<<i<<" "<<mlist[i]<<endl;
  for (int in=1;in<=N;in++) Mmax[in]=Mmax[in-1]+mlist[sps.size()-in];
  // delete [] mlist; //finish up and delete sorted list of spins
  //  for (int in=0;in<=N;in++) cout<<in<<" ["<<Mmin[in]<<", "<<Mmax[in]<<"] "<<endl;

 int *Tmax= new int [N+1];
  Tmax[0]=0; //no particles no spin
  //I have an array sps[i].Jz at my disposal
  // int *mlist=new int [sps.size()];
  for (int i=0;i<sps.size();i++) mlist[i]=sps[i].Tz; //copy list of M
  //QuickSort(sps.size(),mlist); //we sorted the list, print it
  std::sort(mlist, mlist + sps.size());
  //on one level mlist will be -j -j+1 ... j 
  //   for (int i=0;i<sps.size();i++) cout<<i<<" "<<mlist[i]<<endl;
  for (int in=1;in<=N;in++) Tmax[in]=Tmax[in-1]+mlist[sps.size()-in];
  
  /* Sorting for abelian number of quanta */
   int *Wmax= new int [N+1];
   Wmax[0]=0;
    for (int i=0;i<sps.size();i++) mlist[i]=2*sps[i].N+sps[i].L; //create list with quanta
      //QuickSort(sps.size(),mlist); //we sorted the list, print it
      std::sort(mlist, mlist + sps.size());
      for (int in=1;in<=N;in++) Wmax[in]=Wmax[in-1]+mlist[sps.size()-in];
  //  for (int i=0;i<sps.size();i++) cout<<i<<" "<<mlist[i]<<endl;
  //return 1;
  delete [] mlist; //finish up and delete sorted list of spins
 //   for (int in=0;in<=N;in++) cout<<in<<" ["<<", "<<Mmax[in]<<"] "<<endl;


  spsint *currentstate=new spsint [N]; 
 
  //fermions
  if (type==1) {
    //make an initial distribution
    if (sps.size()<N) {cerr<<"fermions do not fit in this space"<<endl; return 0;}
 //this is my current state 
  int *MSUM=new int [N];
  int *TSUM=new int [N];
  int *WSUM=new int [N];
    n=0;
    JZ=0;TZ=0;LL=true; WW=0;
    for (int i=0;i<N;i++) {
      currentstate[i]=i;
      JZ+=sps[i].Jz;
      TZ+=sps[i].Tz;
      LL^=sps[i].P;
      WW+=(sps[i].N<<1)+sps[i].L; //get all quantum numbers initialized
      }
    
  //   cout<<"initial state is done "<<WW<<endl;
    //[xxxx-------] initial distribution
    
    //do an increment
    bool OKflag=true;
    
    while (OKflag) {
      //control conservation laws
      if ((JZ==Jz)&&(TZ==Tz)&&(LL^P)&&(WW<=HW))
     {
       //here the wanted state is ready
       //Numbering starts from 0
       //  cout<<n<<" JZ="<<JZ<<" TZ="<<TZ<<" L="<<LL<<" : ";
       //for (int i=0;i<N;i++) cout<<int(currentstate[i])<<" ";
       //cout<<endl;
       outf.write((char*)(currentstate), N*sizeof(spsint));
	 n++;
     //cerr<<n<<": "<<mbs_id(N, currentstate, sps.size())<<" : "<<int(currentstate[0])<<" "<<int(currentstate[1])<<" "<<sps.size()<<" "<<fmbs_id(N, currentstate, sps.size())<<endl;
#ifdef CSTATUS
         if (S_MBS()) {S_MBS.stream()<<" n=";
         av::PrintInt(S_MBS.stream(),n);
         S_MBS.stream()<<" "<<tav::FractionalConfigurationPosition(sps.size(),N,currentstate)*100.<<"%";
         }
        #endif
	// if ((n%1000000)==0) cout<<" "<<n<<" "<<int(currentstate[0])<<"/"<<sps.size()<<endl;
    //n gives number of states
}
   

      OKflag=false;
      //prepare MSUM;
      MSUM[0]=0;
    for (int in=0;in<N-1;in++) MSUM[in+1]=MSUM[in]+sps[currentstate[in]].Jz;
    TSUM[0]=0;
    for (int in=0;in<N-1;in++) TSUM[in+1]=TSUM[in]+sps[currentstate[in]].Tz;
    WSUM[0]=0;
    for (int in=0;in<N-1;in++) WSUM[in+1]=WSUM[in]+2*sps[currentstate[in]].N+sps[currentstate[in]].L; //2n+l
    //for (int in=0;in<N;in++) cout<<"-------------- "<<MSUM[in]<<endl;
    for (int i=N-1;i>=0;i--) {
  //go in reversed order
      /* I have N-i particles moving starting from I 
      [x-xx-x-Xxx-x-x]
       0 1    i     N-1                 
      say MSUM[i]=sps[0].Jz+sps[1].Jz+..sps[i-1].Jz
      it only worth to proceed if Mmin[N-i]<=Jz-MSUM[i]<=Mmax[N-i]
      */
      
      if (abs(Jz-MSUM[i])>Mmax[N-i]) {continue;}
      if (abs(Tz-TSUM[i])>Tmax[N-i]) {continue;}
     // cerr<<HW<<" "<<WSUM[i]<<" "<<Wmax[N-1]<<endl;
    //if (abs(HW-WSUM[i])>Wmax[N-i]) {continue;}
      // if (abs(Jz+MSUM[i])>Mmax[N-i]) continue;
      
    
      //  if (currentstate[i]+1==sps.size()) continue; //we are at end [x-xxx----X]
      // if ((i!=(N-1))&&((currentstate[i]+1)==currentstate[i+1])) continue;
      // we need to move i-th particle and dump all others after it
      // [x--xX---xxx]->[x--x-Xxxx---]
      //number of spaces I have after it is sps.size()-1-currentstate[i]
      //number of particle including X is N-i
      if ((sps.size()-1-currentstate[i])<(N-i)) continue; //we do not fit;
      //there is a trick here we need to change current i last
      for (int j=N-1;j>=i;j--) {
	//
	//remove all quantum numbers
	JZ-=sps[currentstate[j]].Jz;
	TZ-=sps[currentstate[j]].Tz;
	LL^=sps[currentstate[j]].P;
    WW-=2*sps[currentstate[j]].N+sps[currentstate[j]].L;
	currentstate[j]=currentstate[i]+j+1-i; //move particle
	//put quantum numbers back
	JZ+=sps[currentstate[j]].Jz;
	TZ+=sps[currentstate[j]].Tz;
    WW+=2*sps[currentstate[j]].N+sps[currentstate[j]].L;
	LL^=sps[currentstate[j]].P;

}
      // n++;
//cerr<<WW<<endl;

      OKflag=true;
      break;
    }
    }; 
    delete [] MSUM;
    delete [] TSUM;
    delete [] WSUM;
  }
 //bosons
  if (type==0) {
std::cerr<<"Boson version is not implemented "<<std::endl;

  } //---------------END BOSONS-------------------
  delete [] currentstate;
  //now I need to overwrite n in file
  outf.seekp(0); //go to the beguinning of the file
 outf.write((char*)(&N),sizeof(int));
  outf.write((char*)(&n),sizeof(uint_nbasis)); //overide new n and N
outf.close();
  return n;
}
    /* We make many-body states with rejection by the number of qunta
     */
    uint_nbasis MakeMBS2EQ(
                         const char *filename, ///<output file
                         int N, int Jz, int Tz, int P,
                         int HW, ///< number of oscillator quanta
                         SingleParticleStates &sps ///<single-particle states
    )
    {
        //process status variable if defined
#ifdef CSTATUS
        av::status S_MBS("MBS");
#ifdef FSTATUS
        S_MBS.SetFile();
#endif
#endif
        //-----------end status
        uint_nbasis n=0;
        int type=sps[0].J%2;
        ofstream outf;
        outf.open(filename, ios::binary);
        outf.write((char*)(&N),sizeof(int));
        outf.write((char*)(&n),sizeof(uint_nbasis)); //write zero for now
        
        spin TZ;
        spin JZ;
        bool LL; //abelian spin and parity variables
        int WW; //number of quanta
        cerr<<"Number of quanta "<<HW<<endl;
        //create list Mmax
        //I need a table given that ni particles what maximum M and T they could have
        int *Mmax= new int [N+1];
        Mmax[0]=0; //no particles no spin
        //I have an array sps[i].Jz at my disposal
        int *mlist=new int [sps.size()];
        for (int i=0;i<sps.size();i++) mlist[i]=sps[i].Jz; //copy list of M
        //QuickSort(sps.size(),mlist); //we sorted the list, print it
        std::sort(mlist, mlist + sps.size());
        //on one level mlist will be -j -j+1 ... j
        //   for (int i=0;i<sps.size();i++) cout<<i<<" "<<mlist[i]<<endl;
        for (int in=1;in<=N;in++) Mmax[in]=Mmax[in-1]+mlist[sps.size()-in];
        // delete [] mlist; //finish up and delete sorted list of spins
        //  for (int in=0;in<=N;in++) cout<<in<<" ["<<Mmin[in]<<", "<<Mmax[in]<<"] "<<endl;
        
        int *Tmax= new int [N+1];
        Tmax[0]=0; //no particles no spin
        //I have an array sps[i].Jz at my disposal
        // int *mlist=new int [sps.size()];
        for (int i=0;i<sps.size();i++) mlist[i]=sps[i].Tz; //copy list of M
        //QuickSort(sps.size(),mlist); //we sorted the list, print it
        std::sort(mlist, mlist + sps.size());
        //on one level mlist will be -j -j+1 ... j
        //   for (int i=0;i<sps.size();i++) cout<<i<<" "<<mlist[i]<<endl;
        for (int in=1;in<=N;in++) Tmax[in]=Tmax[in-1]+mlist[sps.size()-in];
        
        /* Sorting for abelian number of quanta */
        int *Wmax= new int [N+1];
        Wmax[0]=0;
        for (int i=0;i<sps.size();i++) mlist[i]=2*sps[i].N+sps[i].L; //create list with quanta
        //QuickSort(sps.size(),mlist); //we sorted the list, print it
        std::sort(mlist, mlist + sps.size());
        for (int in=1;in<=N;in++) Wmax[in]=Wmax[in-1]+mlist[sps.size()-in];
        //  for (int i=0;i<sps.size();i++) cout<<i<<" "<<mlist[i]<<endl;
        //return 1;
        delete [] mlist; //finish up and delete sorted list of spins
        //   for (int in=0;in<=N;in++) cout<<in<<" ["<<", "<<Mmax[in]<<"] "<<endl;
        
        
        spsint *currentstate=new spsint [N];
        
        //fermions
        if (type==1) {
            //make an initial distribution
            if (sps.size()<N) {cerr<<"fermions do not fit in this space"<<endl; return 0;}
            //this is my current state
            int *MSUM=new int [N];
            int *TSUM=new int [N];
            int *WSUM=new int [N];
            n=0;
            JZ=0;TZ=0;LL=true; WW=0;
            for (int i=0;i<N;i++) {
                currentstate[i]=i;
                JZ+=sps[i].Jz;
                TZ+=sps[i].Tz;
                LL^=sps[i].P;
                WW+=(sps[i].N<<1)+sps[i].L; //get all quantum numbers initialized
            }
            
            //   cout<<"initial state is done "<<WW<<endl;
            //[xxxx-------] initial distribution
            
            //do an increment
            bool OKflag=true;
            
            while (OKflag) {
                //control conservation laws
                if ((JZ==Jz)&&(TZ==Tz)&&(LL^P)&&(WW==HW))
                {
                    //here the wanted state is ready
                    //Numbering starts from 0
                    //  cout<<n<<" JZ="<<JZ<<" TZ="<<TZ<<" L="<<LL<<" : ";
                    //for (int i=0;i<N;i++) cout<<int(currentstate[i])<<" ";
                    //cout<<endl;
                    outf.write((char*)(currentstate), N*sizeof(spsint));
                    n++;
                    //cerr<<n<<": "<<mbs_id(N, currentstate, sps.size())<<" : "<<int(currentstate[0])<<" "<<int(currentstate[1])<<" "<<sps.size()<<" "<<fmbs_id(N, currentstate, sps.size())<<endl;
#ifdef CSTATUS
                    if (S_MBS()) {S_MBS.stream()<<" n=";
                        av::PrintInt(S_MBS.stream(),n);
                        S_MBS.stream()<<" "<<tav::FractionalConfigurationPosition(sps.size(),N,currentstate)*100.<<"%";
                    }
#endif
                    // if ((n%1000000)==0) cout<<" "<<n<<" "<<int(currentstate[0])<<"/"<<sps.size()<<endl;
                    //n gives number of states
                }
                
                
                OKflag=false;
                //prepare MSUM;
                MSUM[0]=0;
                for (int in=0;in<N-1;in++) MSUM[in+1]=MSUM[in]+sps[currentstate[in]].Jz;
                TSUM[0]=0;
                for (int in=0;in<N-1;in++) TSUM[in+1]=TSUM[in]+sps[currentstate[in]].Tz;
                WSUM[0]=0;
                for (int in=0;in<N-1;in++) WSUM[in+1]=WSUM[in]+2*sps[currentstate[in]].N+sps[currentstate[in]].L; //2n+l
                //for (int in=0;in<N;in++) cout<<"-------------- "<<MSUM[in]<<endl;
                for (int i=N-1;i>=0;i--) {
                    //go in reversed order
                    /* I have N-i particles moving starting from I
                     [x-xx-x-Xxx-x-x]
                     0 1    i     N-1
                     say MSUM[i]=sps[0].Jz+sps[1].Jz+..sps[i-1].Jz
                     it only worth to proceed if Mmin[N-i]<=Jz-MSUM[i]<=Mmax[N-i]
                     */
                    
                    if (abs(Jz-MSUM[i])>Mmax[N-i]) {continue;}
                    if (abs(Tz-TSUM[i])>Tmax[N-i]) {continue;}
                    // cerr<<HW<<" "<<WSUM[i]<<" "<<Wmax[N-1]<<endl;
                    //if (abs(HW-WSUM[i])>Wmax[N-i]) {continue;}
                    // if (abs(Jz+MSUM[i])>Mmax[N-i]) continue;
                    
                    
                    //  if (currentstate[i]+1==sps.size()) continue; //we are at end [x-xxx----X]
                    // if ((i!=(N-1))&&((currentstate[i]+1)==currentstate[i+1])) continue;
                    // we need to move i-th particle and dump all others after it
                    // [x--xX---xxx]->[x--x-Xxxx---]
                    //number of spaces I have after it is sps.size()-1-currentstate[i]
                    //number of particle including X is N-i
                    if ((sps.size()-1-currentstate[i])<(N-i)) continue; //we do not fit;
                    //there is a trick here we need to change current i last
                    for (int j=N-1;j>=i;j--) {
                        //
                        //remove all quantum numbers
                        JZ-=sps[currentstate[j]].Jz;
                        TZ-=sps[currentstate[j]].Tz;
                        LL^=sps[currentstate[j]].P;
                        WW-=2*sps[currentstate[j]].N+sps[currentstate[j]].L;
                        currentstate[j]=currentstate[i]+j+1-i; //move particle
                        //put quantum numbers back
                        JZ+=sps[currentstate[j]].Jz;
                        TZ+=sps[currentstate[j]].Tz;
                        WW+=2*sps[currentstate[j]].N+sps[currentstate[j]].L;
                        LL^=sps[currentstate[j]].P;
                        
                    }
                    // n++;
                    //cerr<<WW<<endl;
                    
                    OKflag=true;
                    break;
                }
            }; 
            delete [] MSUM;
            delete [] TSUM;
            delete [] WSUM;
        }
        //bosons
        if (type==0) {
            std::cerr<<"Boson version is not implemented "<<std::endl;
            
        } //---------------END BOSONS-------------------
        delete [] currentstate;
        //now I need to overwrite n in file
        outf.seekp(0); //go to the beguinning of the file
        outf.write((char*)(&N),sizeof(int));
        outf.write((char*)(&n),sizeof(uint_nbasis)); //overide new n and N
        outf.close();
        return n;
    }
///Function to test many-body states, saved in file
int TestMBS(char *filename) {
  ifstream outf;
  outf.open(filename, ios::binary);
if (!outf)   {               // if the file does not exist
  FatalError("File "<<filename<<" is not found"<<endl)
return 0; 
}
 int N;
 uint_nbasis n;
  outf.read((char*)(&N),sizeof(int));
  outf.read((char*)(&n),sizeof(uint_nbasis));
  //preparing space for z
  cout<<"Number of particls "<<N<<endl;
  cout<<"Dimension of basis "<<n<<endl;
  spsint *z=new spsint [N];
  spsint *zp=new spsint [N];
  outf.read((char*)(zp), N*sizeof(spsint));
   if (outf.eof()) cerr<<"We reached EOF on line nn=0"<<endl;
   if (outf.fail()) cerr<<"Reading failed on line nn=0"<<endl;
   //----print first line
//   DBG(for (int i=0;i<N;i++) cerr<<int(zp[i])<<" "; cerr<<endl;)
   //-----print end
   int ptest=-3;
   bool ford=true;
   bool norder=true;
  int minz;
   int maxz;
   minz=maxz=zp[0];
   for (int in=1;in<N;in++) { 
     if (zp[in]<minz) minz=zp[in]; 
     if (zp[in]>maxz) maxz=zp[in];}



  for (uint_nbasis nn=1;nn<n;nn++) 
  {
   outf.read((char*)(z), N*sizeof(spsint));
   if (outf.eof()) cerr<<"We reached EOF on line nn="<<nn<<endl;
   if (outf.fail()) cerr<<"Reading failed on line nn="<<nn<<endl;
// -----print
//DBG(for (int i=0;i<N;i++) cerr<<int(z[i])<<" "; cerr<<endl;)
// --print end
 for (int in=0;in<N;in++) { 
     if (z[in]<minz) minz=z[in];
     if (z[in]>maxz) maxz=z[in];}
  for (int in=1;in<N;in++) { 
    if (z[in]<=z[in-1]) {
        std::cerr << "bad bosons at nn="<<nn   <<std::endl<< __FILE__ << ":" << __LINE__ << std::endl;
      norder=false;
    }
  }

 //cerr<<"compare "<<compare(N,zp,z)<<endl;
   if (nn==1)   ptest=compare(N,zp,z);
   else {
     if(ptest!=compare(N,zp,z)) {
       ford=false; 
       std::cerr << "no order at nn="<<nn   <<std::endl<< __FILE__ << ":" << __LINE__ << std::endl;
	   }
   }
   //need swap
   spsint *zz;
   zz=zp; zp=z; z=zz;
   }
  delete [] z;
  delete [] zp;
  if (!ford) cerr<<"States are NOT ordered"<<endl;
 else cerr<<"States are ordered"<<endl;
  if (!norder) cerr<<"Particle order is not present THIS IS FATAL (or bosons)"<<endl;
 else cerr<<"particles are ordered"<<endl;
  cerr<<"SPS are such that N is in ["<<minz<<", "<<maxz<<"]"<<endl;
  outf.close();
  return 1;
}

///Print Basis
uint_nbasis Print (int N, int Jz, int Tz, int P, SingleParticleStates &sps)
{
 int type=sps[0].J%2;
  spin TZ;
  spin JZ;
  bool LL; //abelian spin and parity variables

  uint_nbasis n;
  spsint *currentstate=new spsint [N]; 
   n=0;
    JZ=0;TZ=0;LL=true; //start with opposite parity
  //this is my current state 

  //fermions
  if (type==1) {
    //make an initial distribution
    if (sps.size()<N) {cerr<<"fermions do not fit in this space"<<endl; return 0;}
    for (int i=0;i<N;i++) {
      currentstate[i]=i;
      JZ+=sps[i].Jz;
      TZ+=sps[i].Tz;
      LL^=sps[i].P;}
    
    // cout<<"initial state is done"<<endl;
    //[xxxx-------] initial distribution
    
    //do an increment
    bool OKflag=true;
    
    while (OKflag) {
      //control conservation laws
      if ((JZ==Jz)&&(TZ==Tz)&&(LL^P)) 
     {
       //here the wanted state is ready
       //Numbering starts from 0
      cout<<n<<" JZ="<<JZ<<" TZ="<<TZ<<" L="<<LL<<" : ";
    for (int i=0;i<N;i++) cout<<int(currentstate[i])<<" ";
    cout<<endl;
    n++;
    //n gives number of states
}
   

      OKflag=false;
    for (int i=N-1;i>=0;i--) { 
      //go in reversed order
      //  if (currentstate[i]+1==sps.size()) continue; //we are at end [x-xxx----X]
      // if ((i!=(N-1))&&((currentstate[i]+1)==currentstate[i+1])) continue;
      // we need to move i-th particle and dump all others after it
      // [x--xX---xxx]->[x--x-Xxxx---]
      //number of spaces I have after it is sps.size()-1-currentstate[i]
      //number of particle including X is N-i
      if ((sps.size()-1-currentstate[i])<(N-i)) continue; //we do not fit;
      //there is a trick here we need to change current i last
      for (int j=N-1;j>=i;j--) {
	//remove all quantum numbers
	JZ-=sps[currentstate[j]].Jz;
	TZ-=sps[currentstate[j]].Tz;
	LL^=sps[currentstate[j]].P;
	currentstate[j]=currentstate[i]+j+1-i; //move particle
	//put quantum numbers back
	JZ+=sps[currentstate[j]].Jz;
	TZ+=sps[currentstate[j]].Tz;
	LL^=sps[currentstate[j]].P;

}
      // n++;


      OKflag=true;
      break;
    }
    }; 

  }
 //bosons
  if (type==0) {
    //make an initial distribution
    for (int i=0;i<N;i++) {
      currentstate[i]=0;
      JZ+=sps[0].Jz;
      TZ+=sps[0].Tz;
      LL^=sps[0].P;
    }  
    // cout<<"initial state is done"<<endl;
    //[xxxx-------] initial distribution
    
    //do an increment
    bool OKflag=true;
    
    while (OKflag) {

     //control conservation laws
      if ((JZ==Jz)&&(TZ==Tz)&&(LL^P)) 
     {
       //here the wanted state is ready
       //Numbering starts from 0
      cout<<n<<" JZ="<<JZ<<" TZ="<<TZ<<" L="<<LL<<" : ";
    for (int i=0;i<N;i++) cout<<int(currentstate[i])<<" ";
    cout<<endl;
    n++;
    //n gives number of states
}
   
 
      OKflag=false;
    for (int i=N-1;i>=0;i--) { //go in reversed order
      //  if (currentstate[i]+1==sps.size()) continue; //we are at end [x-xxx----X]
      // if ((i!=(N-1))&&((currentstate[i]+1)==currentstate[i+1])) continue;
      // we need to move i-th particle and dump all others after it
      // [x--xX---xxx]->[x--x-Xxxx---]
      //number of spaces I have after it is sps.size()-1-currentstate[i]
      //number of particle including X is N-i
      if (currentstate[i]==(sps.size()-1)) continue; //we are at end;
      //there is a trick here we need to change current i last
      for (int j=N-1;j>=i;j--)  {
	JZ-=sps[currentstate[j]].Jz;
	TZ-=sps[currentstate[j]].Tz;
	LL^=sps[currentstate[j]].P;
	 //move particle
	currentstate[j]=currentstate[i]+1;
	//put quantum numbers back
	JZ+=sps[currentstate[j]].Jz;
	TZ+=sps[currentstate[j]].Tz;
	LL^=sps[currentstate[j]].P;


      } 
      OKflag=true;
      break;
    }
    }; 

  }
  delete [] currentstate;
return sps.size();
}
    
    uint_nbasis MakePartitionMBS(char *filename,//Temp file that will save mbs.
                                 int N, int Jz, int Tz, int P,//Quantum Numbers for partition
                                 SM::SingleParticleStates &sps,//sps allowed to be used.
                                 spsint* partition,//The partition
                                 int parts//size of partition
    )
    {
        uint_nbasis n=0;
        //int type=sps[0].J%2;
        ofstream outf;
        outf.open(filename, ios::binary);
        outf.write((char*)(&N),sizeof(int));
        outf.write((char*)(&n),sizeof(uint_nbasis)); //write zero for now
        
        spin TZ;
        spin JZ;
        bool LL; //abelian spin and parity variables
        
        std::vector<spsint> testpart(parts,0);
        spsint *currentstate=new spsint [N];
        //this is my current state
        
        //make an initial distribution
        if (sps.size()<N) {cerr<<"fermions do not fit in this space"<<endl; return 0;}
        n=0;
        JZ=0;TZ=0;LL=true;
        for (int i=0;i<N;i++)
        {
            currentstate[i]=i;
            JZ+=sps[i].Jz;
            TZ+=sps[i].Tz;
            LL^=sps[i].P;
        }
        
        // cout<<"initial state is done"<<endl;
        //[xxxx-------] initial distribution
        
        //do an increment
        bool OKflag=true;
        
        while (OKflag)
        {
            //control conservation laws
            if ((JZ==Jz)&&(TZ==Tz)&&(LL^P))
            {
                //here the wanted state is ready
                //Numbering starts from 0
                //  cout<<n<<" JZ="<<JZ<<" TZ="<<TZ<<" L="<<LL<<" : ";
                //for (int i=0;i<N;i++) cout<<int(currentstate[i])<<" ";
                //cout<<endl;
                int indexing=0;
                for (int i=0;i<N;i++)
                {
                    //                    indexing=sps[currentstate[i]].N*2+sps[currentstate[i]].L;
                    //                if (sps[currentstate[i]].Tz<0) indexing+=parts/2;
                    indexing=sps[currentstate[i]].level;
                    testpart[indexing]++;
                }
                bool temptestflag=true;
                for (int i=0;i<parts;i++)
                {
                    if (testpart[i]!=partition[i]) {temptestflag=false; break;}
                }
                if (temptestflag)
                {
                    outf.write((char*)(currentstate), N*sizeof(spsint));
                    n++;
                }
                //            else
                //            {
                //                cout<<"REJECTED: ";
                //                for (int i=0;i<testpart.size();i++) cout<<int(testpart[i])<<" ";
                //                cout<<endl;
                //            }
                std::fill(testpart.begin(), testpart.end(), 0);
                // if ((n%1000)==0) cout<<n<<endl;
                //n gives number of states
            }
            
            
            OKflag=false;
            for (int i=N-1;i>=0;i--)
            {
                //go in reversed order
                //  if (currentstate[i]+1==sps.size()) continue; //we are at end [x-xxx----X]
                // if ((i!=(N-1))&&((currentstate[i]+1)==currentstate[i+1])) continue;
                // we need to move i-th particle and dump all others after it
                // [x--xX---xxx]->[x--x-Xxxx---]
                //number of spaces I have after it is sps.size()-1-currentstate[i]
                //number of particle including X is N-i
                if ((sps.size()-1-currentstate[i])<(N-i)) continue; //we do not fit;
                //there is a trick here we need to change current i last
                for (int j=N-1;j>=i;j--) {
                    //remove all quantum numbers
                    JZ-=sps[currentstate[j]].Jz;
                    TZ-=sps[currentstate[j]].Tz;
                    LL^=sps[currentstate[j]].P;
                    currentstate[j]=currentstate[i]+j+1-i; //move particle
                    //put quantum numbers back
                    JZ+=sps[currentstate[j]].Jz;
                    TZ+=sps[currentstate[j]].Tz;
                    LL^=sps[currentstate[j]].P;
                    
                }
                // n++;
                
                
                OKflag=true;
                break;
            }
        };
        
        
        delete [] currentstate;
        //now I need to overwrite n in file
        outf.seekp(0); //go to the beguinning of the file
        outf.write((char*)(&N),sizeof(int));
        outf.write((char*)(&n),sizeof(uint_nbasis)); //overide new n and N
        outf.close();
        return n;
    }
    uint_nbasis MakePartitionPNMBS(const char *filename,//Temp file that will save mbs.
                                 int N, int Jz, int Tz, int P,//Quantum Numbers for partition
                                 SM::SingleParticleStates &sps,//sps allowed to be used.
                                 spsint* partition,//The partition
                                 int parts//size of partition
    )
    {
        uint_nbasis n=0;
        
        int maxnLev = 0;
        for (int x = 0; x < sps.size(); ++x)
            if (sps[x].Tz > 0 && sps[x].level > maxnLev) maxnLev = sps[x].level;

        maxnLev++;
        //int type=sps[0].J%2;
        ofstream outf;
        outf.open(filename, ios::binary);
        outf.write((char*)(&N),sizeof(int));
        outf.write((char*)(&n),sizeof(uint_nbasis)); //write zero for now
        
        spin TZ;
        spin JZ;
        bool LL; //abelian spin and parity variables
        
        std::vector<spsint> testpart(parts,0);
        spsint *currentstate=new spsint [N];
        //this is my current state
        
        //make an initial distribution
        if (sps.size()<N) {cerr<<"fermions do not fit in this space"<<endl; return 0;}
        n=0;
        JZ=0;TZ=0;LL=true;
        for (int i=0;i<N;i++)
        {
            currentstate[i]=i;
            JZ+=sps[i].Jz;
            TZ+=sps[i].Tz;
            LL^=sps[i].P;
        }
        
        // cout<<"initial state is done"<<endl;
        //[xxxx-------] initial distribution
        
        //do an increment
        bool OKflag=true;
        
        while (OKflag)
        {
            //control conservation laws
            if ((JZ==Jz)&&(TZ==Tz)&&(LL^P))
            {
                //here the wanted state is ready
                //Numbering starts from 0
                //  cout<<n<<" JZ="<<JZ<<" TZ="<<TZ<<" L="<<LL<<" : ";
                //for (int i=0;i<N;i++) cout<<int(currentstate[i])<<" ";
                //cout<<endl;
                int indexing=0;
                for (int i=0;i<N;i++)
                {
                    //                    indexing=sps[currentstate[i]].N*2+sps[currentstate[i]].L;
                    //                if (sps[currentstate[i]].Tz<0) indexing+=parts/2;
                    indexing=sps[currentstate[i]].level;
                    if (sps[currentstate[i]].Tz < 0) indexing += maxnLev;
                    testpart[indexing]++;
                }
                bool temptestflag=true;
                for (int i=0;i<parts;i++)
                {
                    if (testpart[i]!=partition[i]) {temptestflag=false; break;}
                }
                if (temptestflag)
                {
                    outf.write((char*)(currentstate), N*sizeof(spsint));
                    n++;
                }
                //            else
                //            {
                //                cout<<"REJECTED: ";
                //                for (int i=0;i<testpart.size();i++) cout<<int(testpart[i])<<" ";
                //                cout<<endl;
                //            }
                std::fill(testpart.begin(), testpart.end(), 0);
                // if ((n%1000)==0) cout<<n<<endl;
                //n gives number of states
            }
            
            
            OKflag=false;
            for (int i=N-1;i>=0;i--)
            {
                //go in reversed order
                //  if (currentstate[i]+1==sps.size()) continue; //we are at end [x-xxx----X]
                // if ((i!=(N-1))&&((currentstate[i]+1)==currentstate[i+1])) continue;
                // we need to move i-th particle and dump all others after it
                // [x--xX---xxx]->[x--x-Xxxx---]
                //number of spaces I have after it is sps.size()-1-currentstate[i]
                //number of particle including X is N-i
                if ((sps.size()-1-currentstate[i])<(N-i)) continue; //we do not fit;
                //there is a trick here we need to change current i last
                for (int j=N-1;j>=i;j--) {
                    //remove all quantum numbers
                    JZ-=sps[currentstate[j]].Jz;
                    TZ-=sps[currentstate[j]].Tz;
                    LL^=sps[currentstate[j]].P;
                    currentstate[j]=currentstate[i]+j+1-i; //move particle
                    //put quantum numbers back
                    JZ+=sps[currentstate[j]].Jz;
                    TZ+=sps[currentstate[j]].Tz;
                    LL^=sps[currentstate[j]].P;
                    
                }
                // n++;
                
                
                OKflag=true;
                break;
            }
        };
        
        
        delete [] currentstate;
        //now I need to overwrite n in file
        outf.seekp(0); //go to the beguinning of the file
        outf.write((char*)(&N),sizeof(int));
        outf.write((char*)(&n),sizeof(uint_nbasis)); //overide new n and N
        outf.close();
        return n;
    }

//    //This converts the many body states from sps_temp to sps. Result is stored back in the same many
//    //body states object that was given to the function. Assume no sorting is neccessary. (this is due to the
//    //fact that this will be used when some shells are just omitted so ordering is preserved).
//#ifndef ConvertMBSFUNC
//#define ConvertMBSFUNC
//    void ConvertMBS(SM::Many_Body_States &st, SM::SingleParticleStates &sps,SM::SingleParticleStates &sps_temp)
//    {
//        vector<spsint> key(sps_temp.size());
//        spsint tempkey;
//        for (int j=0;j<sps_temp.size();j++)
//        {
//            for (spsint i=0;i<sps.size();i++)
//            {
//                if (sps_temp[j].N==sps[i].N && sps_temp[j].L==sps[i].L && sps_temp[j].J==sps[i].J && sps_temp[j].Jz==sps[i].Jz && sps_temp[j].Tz==sps[i].Tz) {tempkey=i;break;}
//            }
//            //            key[i]=SpsIndex(sps_temp,sps,i);
//            key[j]=tempkey;
//        }
//        
//        for (uint_nbasis i=0;i<st.n;i++)
//        {
//            for (int j=0;j<st.N;j++)
//                st.z[i][j]=key[st.z[i][j]];
//        }
//    }
//#endif
//---------------------------------------------------------------------
} //end namespace sm

#endif //__MAKEMBS__CXX__

