/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 

\author Alexander Volya  <http://www.volya.net>
\file statusclass.cxx
\brief Ouput status of running programs, OMP-master thread.
\author A. Volya
*/
/*
03/17/09 PrintTime Function is introduced for nice output
*/
#ifndef _STATUSCLASS_CXX_
#define _STATUSCLASS_CXX_
#include <iostream>
#include <fstream>
using std::ofstream;
//#include <sstream>
//using std::ostringstream;
#include <iomanip>
using std::setw;
using std::setfill;
#include <cstring>
//using std::string;
using std::strcpy;
//#include <cstdlib>
//using std::setprecision
#include <ctime>
using std::cerr;
using std::clock;
namespace av{
#define MYTIME_IN_SEC (double(clock())/CLOCKS_PER_SEC)

//********************************************************************
/** \brief print time, output 7+6=13 characters

\todo return does not work, it would be logical to have spearate print functions with stringstreams 
*/
std::ostream&  PrintTime(std::ostream &logfile, double itime=MYTIME_IN_SEC) {
logfile<<std::fixed<<std::setprecision(2);
double mytime=itime;
if (mytime<60.) {logfile<<setw(7)<<mytime<<" sec  "<<std::setprecision(6); return logfile;} 
mytime/=60.;
if (mytime<60.) {logfile<<setw(7)<<mytime<<" min  "<<std::setprecision(6); return logfile;}
mytime/=60.;
if (mytime<24.) {logfile<<setw(7)<<mytime<<" hr   "<<std::setprecision(6); return logfile;}
mytime/=24.;   
logfile<<setw(8)<<mytime<<" days "<<std::setprecision(6); return logfile;
}

std::ostream&  PrintInt(std::ostream &logfile, unsigned long ix) {
logfile<<std::fixed<<std::setprecision(2);
double x=ix;
if (x<1000.) {logfile<<setw(7)<<ix<<"    "<<std::setprecision(6); return logfile;} 
x/=1000.;
if (x<1000.) {logfile<<setw(7)<<x<<" k  "<<std::setprecision(6); return logfile;}
x/=1000.;
if (x<1000.) {logfile<<setw(7)<<x<<" M  "<<std::setprecision(6); return logfile;}
x/=1000.;
if (x<1000.) {logfile<<setw(7)<<x<<" G  "<<std::setprecision(6); return logfile;}
x/=1000.;
if (x<1000.) {logfile<<setw(7)<<x<<" T  "<<std::setprecision(6); return logfile;}
x/=1000.;
if (x<1000.) {logfile<<setw(7)<<x<<" P  "<<std::setprecision(6); return logfile;}
x/=1000.;
if (x<1000.) {logfile<<setw(7)<<x<<" E  "<<std::setprecision(6); return logfile;}
return logfile;
}


class status  {
	private:
	time_t reset_time; ///<time at which object is reated
	time_t last_check_time; ///<time at which last status was checked
	time_t present_time; ///< time for present
	time_t update_interval_time_max; ///<update interval in units of time maximum
	time_t update_interval_time_min; ///<low value, increase
	time_t reset_time_clock; ///<this time is given with clock
	double clock_rate; ///<internal clock rate
	double run_time;  ///<running time in seconds
	unsigned long progress; ///<progress as an ascending integer
	unsigned long progress_step; ///<step in progress
	unsigned long progress_end; ///< ending integer
	unsigned long progress_next_print; ///<progress+progress_step, next time to update
	//unsigned long update_interval_sec; ///< update interval in seconds
	unsigned int update_bitshift_int; ///< for integer iterations
	std::ofstream filestream; ///<Internal file stream if we save to file
    std::ostream *statusstream; ///<pointer to a stream used for output
	char cret; ///<carriage return character
//functions

/*! \brief Internal initialization function
*/
void Initialize (
const double min_update_time, ///< update frequencey
const char* _name ///< [in] name of the status variable
) {
strcpy(name,_name);
reset_time=time(NULL);
last_check_time=reset_time;
present_time=reset_time;
SetRefreshTime(min_update_time);
reset_time_clock=clock();
update_bitshift_int=1; //initial update shift is 0 
progress_step=1;
progress_end=0; //default do not know when we end
}

/* If access to this function happens too often we increase update_bitshift 
in the case of too infrequent access we decrease update_bitshift. 
computes run_time and clock_rate
*/
void BalanceInterval(void) {
present_time=time(NULL);
//std::cerr<<(present_time-last_check_time)<<" "<<update_interval_time_min<<" "<<progress_step<<std::endl;
if ((present_time-last_check_time)<update_interval_time_min) progress_step<<=1; //rescale timing
if ((present_time-last_check_time)>(update_interval_time_max)) if (progress_step>1) progress_step>>=1; //rescale timing back 
 last_check_time=present_time;
//do other things, set rate and present time
 progress=present_time; //set progress based on time
 run_time=(present_time-reset_time);
time_t present_time_clock=clock();
double run_time_clock=double(present_time_clock-reset_time_clock)/double(CLOCKS_PER_SEC);
if (run_time!=0) clock_rate=run_time_clock/run_time; else clock_rate=1.0;  

}

/*! \brief Process status using time
*/
bool ProcessStatus(void) {
unsigned long current_progress=time(NULL);
//if ((progress>>update_bitshift_int)<<update_bitshift_int==progress) std::cerr<<" progres "<<progress<<std::endl;
if (current_progress>=progress_next_print) {
BalanceInterval();
progress_next_print=progress+progress_step; //next output is printed when progress> progress_next_print
return true;
} else return false;
}

/*! \brief Process status using integer
*/
bool ProcessStatus(const long n) {
//if ((progress>>update_bitshift_int)<<update_bitshift_int==progress) std::cerr<<" progres "<<progress<<std::endl;
if (n>=progress_next_print) {
BalanceInterval(); progress=n;
progress_next_print=progress+progress_step; //next output is printed when progress> progress_next_print
return true;
} else return false;
}


public:
char name[512]; //status name
//operator std::ostream&() {return (*statusstream);} ///< This casts class to stream for output
std::ostream& stream(void) {return (*statusstream);}  ///< returns the stream
/**! \brief Constructor
sets last reset time and last check time
*/


void Reset(const char* _name="") {
Initialize(1.0,_name);
progress=present_time;
progress_next_print=progress+progress_step;
if (!filestream.is_open()) {
statusstream=&std::cerr; ///set default output stream
cret='\r'; ///<set default carriage return
}
(*statusstream)<<std::endl;
(*statusstream)<<name<<": ";
(*statusstream)<<"Starting: "<<name;
(*statusstream)<<ctime (&reset_time);
(*statusstream).flush();
//progress=reset_time; //progress is based on time
}

status (const char* _name="")  {
Reset(_name);}

void Reset (long n,const char* _name="")  {
Initialize(1.0,_name);
progress=0; //progress is based on time
progress_next_print=progress+progress_step;
progress_end=n;
if (!filestream.is_open()) { //if file is open continue writing
statusstream=&std::cerr; ///set default output stream
cret='\r'; ///<set default carriage return
}
(*statusstream)<<std::endl;
(*statusstream)<<name<<": ";
(*statusstream)<<"Starting: ";
(*statusstream)<<ctime (&reset_time);
(*statusstream).flush();
}

status (long n,const char* _name="") {Reset(n,_name);}




void SetFile(const char* fname) {
if (filestream.is_open()) filestream.close();
filestream.open(fname);
std::cerr<<"Status saved in the log file: "<<fname<<std::endl;
statusstream=&filestream;
cret='\n';
(*statusstream)<<std::endl;
(*statusstream)<<ctime (&reset_time);
(*statusstream)<<"Starting: "<<name<<std::endl;
}

void SetFile() {
char fname[256];
snprintf(fname,256,"%s_%i.status", name, int(time (NULL)));
SetFile(fname);
}


~status() {
(*statusstream)<<std::endl;
(*statusstream)<<name<<": ";
(*statusstream)<<"Completed: ";
PrintTime( (*statusstream),run_time);
(*statusstream)<<" RATE="<<std::setprecision(2)<<clock_rate<<" ";
(*statusstream)<<ctime (&present_time);
(*statusstream)<<std::setprecision(6); //unset 
(*statusstream).flush();
if (filestream.is_open()) filestream.close();
}


/*! \Show status based time
 
 \return true if status report was done
 */
bool ShowStatus(void) {
bool bret=false;
#pragma omp master
if (ProcessStatus()) {
(*statusstream)<<cret;
(*statusstream).flush();
(*statusstream)<<name<<":";
PrintTime((*statusstream),run_time); 
(*statusstream)<<" RATE="<<std::setprecision(2)<<clock_rate;
  (*statusstream)<<std::setprecision(6);
bret=true;
}
return bret;
}
    
/*! \Show status based on integer
 
 \return true if status report was done
 */
bool ShowProgress(long i) {
bool bret=false;
#pragma omp master
if (ProcessStatus(i)) {
	(*statusstream)<<cret;
	(*statusstream).flush();
(*statusstream)<<name<<":";
	//print percentage done
    (*statusstream)<<std::fixed<<std::setprecision(2)<<std::setw(7)<<progress*100./progress_end<<"% "<<std::setprecision(6);
	//print current time 
	PrintTime( (*statusstream),run_time); 
	 (*statusstream)<<"ETA: ";
	 double est_time=1.0E5;
	if ((run_time>1)&&(progress>1)) est_time=run_time*(progress_end/double(progress)-1.0);
	PrintTime( (*statusstream),est_time);
	PrintInt((*statusstream),progress); 
	(*statusstream)<<" RATE="<<std::setprecision(2)<<clock_rate;
    (*statusstream)<<std::setprecision(6);
	//std::cerr<<std::endl;
	//std::cerr<<setw(9+13+5+13)<<setfill('\b')<<""<<setfill(' ');
    bret=true;
}
    return bret;
}


bool operator () (long i) {return ShowProgress(i);} ///< Main call function int-based
bool operator () (void) {return ShowStatus();} ///< Show status function time-based

/*! \brief Set pointer of an internal stream 
*/
void SetStream(std::ostream &outstream) {
statusstream=&outstream;
}

void SetRefreshTime(const double min_update_time) {
update_interval_time_min=long (min_update_time); //10 sec for min
update_interval_time_max=(update_interval_time_min<<1); //set initial update at 100 sec for max
}

}; //end class definitions

} //namespace av


#endif


