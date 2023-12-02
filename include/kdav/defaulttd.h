#ifndef  __DEFAULTTD_H
#define __DEFAULTTD_H

//If programme has to be compilates in LINUX/UNIX system operation.
#define OS_LINUX

//If programme has to be compilates in DOS/WINDOWS system operation.
//#define OS_WINDOWS
 
#ifdef OS_LINUX
 using namespace std;
#endif

typedef double ftyp;
#include <cstdlib>
typedef int64_t ntyp;

#endif  //__DEFAULTTD_H





