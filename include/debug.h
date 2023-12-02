/*! \file debug.h 
\brief Debuggind/utility preprocessor functions/settings 
\example debug_example.cpp
\brief Example file how to used debug.h preprocessing
 */

#ifndef __DEBUG_H__
#define __DEBUG_H__
#include <iostream>
#include <cstring> //use for strcpy
#include <limits> //pause
#include <sstream> //convert stream to string in CriticalError
#include <exception>
char debug_cxx_message[1024]; //used in CriticalError
#define BLAME { \
 std::cerr<<""<<__FILE__<<"\n";\
 std::cerr<<"Author: Alexander Volya, http://www.volya.net\n"; \
 std::cerr<<"Version date: "<<__DATE__<<endl;\
 std::cerr<<"Compiler: "<<__VERSION__<<"\n\n"; \
  }
#ifdef DEBUG


#define DBGPrint(x) std::cerr << #x << "=" << x << " in "  << __FILE__ << ":" << __LINE__ << std::endl
/*! \brief Print Debug message */
/*! \brief Run a set of Debug commands  */
#define DBG(commands) commands
#else
/*! 
\brief If \a DBG is defined, print variable x, file name, line number
*/
#define DBGPrint(x)
/*! 
\brief If \a DBG is defined, print variable message, file name, line number
*/
/*! 
\brief If \a DBG is defined, run commands 
*/
#define DBG(commands) 
#endif


#ifdef WARNING

//! \brief If \a WARNING is defined varable, file name and line will be printed
#define WPrint(x) std::cerr << #x << "=" << x << " in "  << __FILE__ << ":" << __LINE__ << std::endl
//! \breif If \a WARNING is defined message, file line and line will be printed
#define WMessage(message) std::cerr << message   <<std::endl<< __FILE__ << ":" << __LINE__ << std::endl
//! \breif If \a WARNING is defined commands will be performed
#define WARNING(commands) commands
#else
//! \brief If \a WARNING is defined varable, file name and line will be printed
#define WPrint(x)
//! \breif If \a WARNING is defined message, file line and line will be printed
#define WMessage(message)
//! \breif If \a WARNING is defined commands will be performed
#define WARNING(commands) 
#endif



#ifdef INFO


#define INFPrint(x) std::cerr << #x << "=" << x << " in "  << __FILE__ << ":" << __LINE__ << std::endl
#define INFMessage(message) std::cerr << message   <<std::endl<< __FILE__ << ":" << __LINE__ << std::endl
#define INF(commands) commands
#else
//! if \a INF is defined, print variable, file name and line.
#define INFPrint(x)
//! if \a INF is defined, print message, file name and line.
#define INFMessage(message)
//! if \a INF is defined, run commands.
#define INF(commands) 
#endif

//These Are Always defined
/*! \def FatalError(message)
\brief Fatal error stop. Print message, file name and line number. 
 */
#define FatalError(message)		 { std::cerr <<"Error: "<< message <<std::endl<< __FILE__ << ":" << __LINE__ << std::endl; throw (1); }
//! \brief Ask user to enter variable \a x from std input 
//#define Ask(x)	  {std::cerr<< __FILE__ << ":" << __LINE__ << std::endl<<"Input "<< #x <<":"; cin>>x;}
#define Ask(x)	  {std::cerr<<"Input "<< #x <<":"; std::cin>>x;}
#define AskUser(y,x)	  {std::cerr<<y; std::cin>>x;}
//#define CriticalError(message){class myexception: public std::exception{virtual const char* what() const throw(){return message;}} myex;throw myex;}
#define CriticalError(message){ \
std::stringstream outs; outs<<message; std::string stmp=outs.str(); std::strcpy(debug_cxx_message,stmp.c_str()); \
class CriticalErrorException: public std::exception{virtual const char* what() const throw(){return debug_cxx_message;}} myex; throw myex; } 

#define CheckInputStream(DataStream_){ \
if (!DataStream_.good()) {        \
if (DataStream_.fail())            \
{class CriticalErrorException: public std::exception{virtual const char* what() const throw() \
{return "Error: Loss of stream integrity while reading \n";}} myex;throw myex;}                \
if (DataStream_.eof())            \
{class CriticalErrorException: public std::exception{virtual const char* what() const throw() \
{return "Error: End-of-file reached while reading \n ";}} myex;throw myex;} \
} \
}
//!\brief pause processing until EOF or \n is encountered
int Pause(void) {
std::cin.ignore();//Ignores previous input, get cin.get working even with new line inputs.
std::cin.get();
 return 1;
}
//!\brief print message and pause processing until EOF or \n is encountered
int Pause(const char *a) {
  std::cerr<<a<<std::endl;
 return Pause();
}

#ifndef CRITICAL___STOP___FLAG
#define CRITICAL___STOP___FLAG
//! \brief This flag determines if critical stop pauses the program (1/0)=yes/no
int critical_stop_flag=1; //this flag determines if we stop or not
#endif
/*! \brief Critical function, print message and pause program 

The program will be paused with \a message printed. The user has following options
- press Enter= \cr to proceed, the program will stop again if next Critical Stop is encountered. 
- press Escape following with Enter to proceed and ignore all future critical stops. This sets #critical_stop_flag=0
- press ctrl+c to terminate. 
 */
#define CriticalStop(message){ if (critical_stop_flag) {std::cerr <<"Critical Stop: ";std::cerr<<message<<'\n';std::cerr<<"Press \n''enter'' to continue, \n''escape+enter'' to ignore future stops, \n''ctrl-c'' to abort "; char ccc___;do{ std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); ccc___ = std::getchar();if (ccc___ == EOF) break;if (ccc___ == '\033') {critical_stop_flag=0; break;}} while (ccc___ != '\n');}}
#endif // __DEBUG_H__ //
