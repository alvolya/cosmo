/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file basic-rw.cxx
 \brief Read and Write functions for the basic fundamental data types
\ingroup gp_rw
*/
#ifndef __BASIC_RW__CXX__
#define __BASIC_RW__CXX__

#include <debug.h> //FatalError depends on this
#include <fstream> //read write require fstream
using std::cout;
using std::endl;
using std::cin;
using std::cerr;
/*! \namespase rw
\brief Binary and text Read/Write operations.
 
 All operations are programmed as Read(variable, file/stream).
 */
namespace rw{

  
    bool TestRWFile(const char *fname) {
    bool ok = static_cast<bool>(std::ofstream(fname).put('a')); // create file
    if(!ok) {
    //std::perror("Error creating file"); return 1;
    return ok;
    }
    std::remove(fname); //delete file
    bool failed = !std::ifstream("fname");
    return failed; //true is good, false is bad
    //if(failed) { std::perror("Error opening deleted file"); return 1;}
        
           }
    
    bool TestFile(const char *fname) {
        
        std::ifstream inf;
        inf.open(fname);
        bool stat=inf.is_open();
        if (stat) inf.close();
       // else {std::cerr<<"Fail to open file "<<fname<<endl;}
        return stat;
    }
    

    int Read(double &T, std::ifstream &inf) {
        //    cerr<<"call read basic\n"<<endl;
        inf.read((char*)(&T),sizeof(double));
        CheckInputStream(inf);
        return 0;
    }
    
    int Write(const double T,std::ofstream &outf) {
       //        cerr<<"call write basic<double> "<<T<<endl;
        outf.write((char*)(&T),sizeof(double)); //write our number;
        if (!outf.good()) CriticalError("Unable to write data stream ");
        return 0;
    } 
    
    int Read(int &T, std::ifstream &inf) {
      //  cerr<<"call read basic\n"<<endl;
        inf.read((char*)(&T),sizeof(int));
        //integrity of read
        CheckInputStream(inf);
        return 0;
    }
    
    int Write(const int T,std::ofstream &outf) {
     //   cerr<<"call write basic <int>"<<T<<endl;
        outf.write((char*)(&T),sizeof(int)); //write our number;
        if (!outf.good()) CriticalError("Unable to write data stream ");
        return 0;
    } 

    int Read(unsigned int &T, std::ifstream &inf) {
        //  cerr<<"call read basic\n"<<endl;
        inf.read((char*)(&T),sizeof(unsigned int));
        //integrity of read
        CheckInputStream(inf);
        return 0;
    }
    
    int Write(const unsigned int T,std::ofstream &outf) {
        //   cerr<<"call write basic <int>"<<T<<endl;
        outf.write((char*)(&T),sizeof(unsigned int)); //write our number;
        if (!outf.good()) CriticalError("Unable to write data stream ");
        return 0;
    } 
    int Read(long &T, std::ifstream &inf) {
        //  cerr<<"call read basic\n"<<endl;
        inf.read((char*)(&T),sizeof(long));
        //integrity of read
        CheckInputStream(inf);
        return 0;
    }
    
    int Write(const long T,std::ofstream &outf) {
        //   cerr<<"call write basic <int>"<<T<<endl;
        outf.write((char*)(&T),sizeof(long)); //write our number;
        if (!outf.good()) CriticalError("Unable to write data stream ");
        return 0;
    } 
    
    int Read(unsigned long &T, std::ifstream &inf) {
        //  cerr<<"call read basic\n"<<endl;
        inf.read((char*)(&T),sizeof(unsigned long));
        //integrity of read
        CheckInputStream(inf);
        return 0;
    }
    
    int Write(const unsigned long T,std::ofstream &outf) {
        //   cerr<<"call write basic <int>"<<T<<endl;
        outf.write((char*)(&T),sizeof(unsigned long)); //write our number;
        if (!outf.good()) CriticalError("Unable to write data stream ");
        return 0;
    } 
    
    int Read(long long &T, std::ifstream &inf) {
        //  cerr<<"call read basic\n"<<endl;
        inf.read((char*)(&T),sizeof(long long));
        //integrity of read
        CheckInputStream(inf);
        return 0;
    }
    
    int Write(const long long T,std::ofstream &outf) {
        //   cerr<<"call write basic <int>"<<T<<endl;
        outf.write((char*)(&T),sizeof(long long)); //write our number;
        if (!outf.good()) CriticalError("Unable to write data stream ");
        return 0;
    } 
    
    int Read(unsigned long long &T, std::ifstream &inf) {
        //  cerr<<"call read basic\n"<<endl;
        inf.read((char*)(&T),sizeof(unsigned long long));
        //integrity of read
        CheckInputStream(inf);
        return 0;
    }
    
    int Write(const unsigned long long T,std::ofstream &outf) {
        //   cerr<<"call write basic <int>"<<T<<endl;
        outf.write((char*)(&T),sizeof(unsigned long long)); //write our number;
        if (!outf.good()) CriticalError("Unable to write data stream ");
        return 0;
    } 

    
    int Read(float &T, std::ifstream &inf) {
        //  cerr<<"call read basic\n"<<endl;
        inf.read((char*)(&T),sizeof(float));
        //integrity of read
        CheckInputStream(inf);
        return 0;
    }
    
    int Write(const float T,std::ofstream &outf) {
        //   cerr<<"call write basic <int>"<<T<<endl;
        outf.write((char*)(&T),sizeof(float)); //write our number;
        if (!outf.good()) CriticalError("Unable to write data stream ");
        return 0;
    } 
}




#endif //__BASIC-RW__CXX__

