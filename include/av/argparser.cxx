/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file argparser.cxx
\author A. Volya 
\ingroup gp_nr
*/
#ifndef __ARGPARSER__CXX__
#define __ARGPARSER__CXX__
#include <debug.h>
#include <cstdlib>
using std::atoi;
#include <string>
#include <vector>
#include <iomanip> // std::put_time
#include <chrono> //c++11 std::chrono::system_clock
//using std::chrono;
using std::string;
namespace av {




/*! \brief Argument class/structure saves help/format/requirements for each argument flag
*/
class argumentflag {
public:
string S_POSIX;
string S_GNU;
string HELP;
int min; //minimum number of times we can have this element
int max; //maximum number of times this element can appear
int id; //id number 0,1,2... that identifies this argument
};
std::ostream& operator << (std::ostream &stream,
                               const argumentflag &T) 
    {
        stream<<"-"<<T.S_POSIX<<"\t--"<<T.S_GNU<<"\t";
        if (T.min>0) stream<<" required";
        if (T.min>1) {
        if (T.max==T.min ) stream<<" parameters= "<<T.min; //<exact number of parameters
        if (T.max<T.min ) stream<<" parameters >"<<T.min; //<limited from below
        if (T.max>T.min) stream<<T.max<<"> parameters >"<<T.min; //<limited from below
        }
        stream<<std::endl;
        stream<<"\t"<<T.HELP;
     return stream;
    };
/*! \brief Main parser function
We have a vector of flags and the corresponding vector of arguments
*/
class argparser  {
std::vector<argumentflag> argx;///< array of arguments
std::vector<string> comline; ///< copy of the command line
std::vector<std::vector<string> > data;
string PRE_POSIX;
string PRE_GNU;
string DESC;
public:
/*Search list of arguments and if not found return the size */
int FindFlag(string sflag ///<Flag POSIX first GNU second
) {
for (int iy=0;iy<argx.size();iy++) if (sflag==argx[iy].S_POSIX) { return argx[iy].id;}
for (int iy=0;iy<argx.size();iy++) if (sflag==argx[iy].S_GNU) {return argx[iy].id;}
return argx.size(); //if not found return size
}
  int SetDescription(string dd) {DESC=dd; return 0;};

/*! \brief Add new argument */
int SetFlag(
string op, ///< Short POSIX format
string og, ///< Long GNU format
string oh, ///< HELP string
int _min=0, ///< Minimum number the argument is encountered, note if minimum is zero and max is not 1 then no flag is also allowed
int _max=-1 ///< Maximum number the argument is encountered (unlimited if negative)
) {
if (FindFlag(op)!=argx.size()) {
int ix=FindFlag(op);
argx[ix].S_POSIX=op; argx[ix].S_GNU=og; argx[ix].HELP=oh;
argx[ix].min=_min; argx[ix].max=_max;
//std::cerr<<"Flag \""<<op<<"\" is already defined \n";
return argx.size();
}
//if (FindFlag(og)!=argx.size()) {std::cerr<<"Flag \""<<og<<"\" is already defined \n"; return argx.size();}
argumentflag v; v.S_POSIX=op; v.S_GNU=og; v.HELP=oh;
v.min=_min; v.max=_max;
//set id for this argument
v.id=argx.size(); //set id for this argument
argx.push_back(v);
return argx.size();
}

int SetDefault(string oh, int _min=0,int _max=-1) {
argx[0].HELP=oh;
argx[0].min=_min;
argx[0].max=_max;
return argx.size();
}

int PrintHelp(){
std::cerr<<DESC<<std::endl;
std::cerr<<"Usage: "<<comline[0]<<" defaults flag [parameters] flag [parameters] ..."<<std::endl;
std::cerr<<"defaults \n\t"<<argx[0].HELP<<std::endl;
for (int i=1;i<argx.size();i++) std::cerr<<argx[i]<<std::endl;
 return 0;
}

/*! \brief Constructor allowed only from the list of arguments*/
argparser(int _argc, char** _argv) {
/* 
Print blame
*/
 std::cerr<<"Program: "<<_argv[0]<<"\n";\
 std::cerr<<"Author: Alexander Volya, http://www.volya.net\n"; \
 std::cerr<<"Version date: "<<__DATE__<<std::endl;\
 std::cerr<<"Compiler: "<<__VERSION__<<"\n\n"; \
/* Read command line arguments*/
for (int i=0;i<_argc;i++) comline.push_back(_argv[i]);
/* Define format for argument flag prefix */
PRE_POSIX="-";
PRE_GNU="--";
/* Set default */
SetFlag("","",""); //set default option id is zero
SetFlag("h","help","Print help",0,1);
}

int CheckArguments() {
//---------------------->First check if Help is set<--------------------------------
if (data[1].size()>0) {PrintHelp(); return 1;}
bool ferrors=false; //this is a trigger that is turned on if there are any errors
//----------------------->here we check default arguments<-----------------------
if ((argx[0].min>0)&&(data[0].size()<argx[0].min)) {
std::cerr<<"Error: Incorrect number of arguments; \nyou must provide at least "<<argx[0].min-1<<" default argument(s)."<<std::endl;
ferrors=true; };
if (data[0].size()>argx[0].max) {
std::cerr<<"Error: Incorrect number of arguments; \nIt is expected that you have not more than "<<argx[0].max-1<<" default argument(s)."<<std::endl;
ferrors=true; };

//check non-default arguments
for (int i=1;i<data.size();i++) {
if (data[i].size()<argx[i].min) {
std::cerr<<"Error: Incorrect number of arguments under flag \"-"<<argx[i].S_POSIX<<"\"; \nyou must provide at least "<<argx[i].min<<" argument(s) with this flag"<<std::endl; ferrors=true; };
if ((argx[i].min==0)&&(argx[i].max!=1)&&(data[i].size()==1)) {
std::cerr<<"Error: at least 1 argument should follow the  \"-"<<argx[i].S_POSIX<<"\" flag"<<std::endl; ferrors=true; };
}
if (ferrors) {
std::cerr<<"Help: "<<std::endl;
PrintHelp();
}
return ferrors;
}

/*! \brief Main function that parses the command line */
int ParseArguments() {
data.resize(argx.size()); //set data size to fit all arguments
int ix=0; //this is an id for an argument variable
//loop over all command line arguments
data[0].push_back(comline[0]); //the program name itself is the default flag 
for(int i=1;i<comline.size();i++) {
//------------>Check if an argument is a command <----------------------------------------
if (
(PRE_POSIX == comline[i].substr (0,PRE_POSIX.size())) || (PRE_GNU == comline[i].substr (0,PRE_GNU.size()))
) {
ix=-1;
for (int iy=0;iy<argx.size();iy++) if (comline[i]==PRE_POSIX+argx[iy].S_POSIX) {ix=argx[iy].id; break;}
for (int iy=0;iy<argx.size();iy++) if (comline[i]==PRE_GNU+argx[iy].S_GNU) {ix=argx[iy].id; break;}
if (ix==-1) std::cerr<<"Warning: unknown flag ignored: "<<comline[i]<<std::endl;
//else {std::cerr<<"using "<<comline[i]<<std::endl;}
//std::cerr<<"i="<<i<<std::endl;
//continue; //skip the flag
} //end prefix check------------------
if (ix==-1) continue; //skip all arguments that follow an anknown flag
//Skip argument if maximum is reached
if (argx[ix].max==data[ix].size()) {
std::cerr<<"Warning: argument \""<<comline[i]<<"\" is treated as default "<<std::endl;
//std::cerr<<ix<<" "<<argx[ix].S_POSIX<<" "<<data[ix].size()<<std::endl;
data[0].push_back(comline[i]);
}
else  data[ix].push_back(comline[i]);
//std::cout<<comline[i]<<" "<<ix<<std::endl;
} //end loop
return CheckArguments();
}

/*!\brief Get an argument */
string operator () (
string flg, ///< flag
int i=0 ///< id
) {
int id=FindFlag(flg);
if (id==argx.size()) FatalError("flag is not in the list of arguments");
if (i>=data[id].size()) FatalError("The "<<i<<"-th  argument with flag -"<<argx[id].S_POSIX<<" is not defined ");
return data[id][i];
}

bool isSet (
string flg, ///< flag
int i=0 ///< id
) {
int id=FindFlag(flg);
if (id==argx.size()) return false;
if (i>=data[id].size()) return false;
return true;
}

/*number of arguments of this flag
*/
int nArgs (
string flg ///< flag
) {
int id=FindFlag(flg);
if (id==argx.size()) return 0; //no arguments of this type allowd
return data[id].size();
}



/*!\brief Return default argument
*/
string operator () (int i) {
int id=0; ///
if (i>=data[id].size()) FatalError("The "<<i<<"-th  default argument is not defined ");
return data[id][i];
}

/*! \brief return a specific argument of a given flag
*/
string operator () (
string flg, ///<flag
string fdefault, ///<default value
int i=0 ///<number of the argument under this flag
) {
int id=FindFlag(flg);
if (id==argx.size()) FatalError("flag is not in the list of arguments");
string rstring=fdefault;
if (i<data[id].size()) rstring=data[id][i];
return rstring;
}

int intarg(string flg, int i=0) {
return atoi((operator()(flg,i)).c_str());
}

float floatarg(string flg, int i=0) {
return atof((operator()(flg,i)).c_str());
}

std::vector<string> operator [] (string flg) {
int id=FindFlag(flg);
if (id==argx.size()) FatalError("flag is not in the list of arguments");
return data[id];
}

//allow access to private data from cout
friend std::ostream& operator << (std::ostream &stream,
                               const argparser &T);

};

std::ostream& operator << (std::ostream &stream,
                               const argparser &T)
    {
//print current time
auto now = std::chrono::system_clock::now();
auto now_c = std::chrono::system_clock::to_time_t(now);
stream << std::put_time(std::localtime(&now_c), "%c") << " ";
//print command line
stream<<"; ";
for (auto i=0;i<T.comline.size();i++)
stream<<T.comline[i]<<" ";
     return stream;
    };
    
} //end namespace av
#endif //__ARGPARSER__CXX__

