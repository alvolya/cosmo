/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! \file inquire.cxx
*/
#ifndef __INQUIRE__CXX__
#define __INQUIRE__CXX__
#include <iostream>
using std::cerr; 
using std::cin;
#include <string>
using std::string;
#include <sstream>
using std::stringstream;
#include <vector> //inquire from list
using std::vector;

namespace tav {
/*! \brief ask user for variable x with question in variable tx, present default option
*/
template <class number, class dnumber>
number& Inquire(
number &x, //<variable
const char *tx, //<question to user
const dnumber default_number //<default answer
){
cerr<<tx<<"["<<default_number<<"]";
string st;
//std::cin.ignore(256, '\n' ); //ignore all \n characters
//cin.clear();// clear if previous operaton was bad
std::getline(std::cin,st); //get line
//separate valuable substring
 // pos = str.find("live");    // position of "live" in str
size_t found=st.find_first_not_of(" ,.:\t");
if (found!=string::npos) st = st.substr(found);
//char tchar=cin.peek(); //ask for char
if (st.empty()) {
//here we do not know if this is default or not
//std::numeric_limits<std::streamsize>::max() limits are not defined
cerr<<"<>";
//std::cin.ignore(256, '\n' ); //ignore all \n characters
std::getline(std::cin,st); //get line
}
//else std::getline(std::cin,st); //get line
  
//long q=cin.tellg(); //this appears to remove '/n' if presnt from the cin buffer Very strage but works!


stringstream ss;
ss<<st;
if (st.empty()) x=default_number;
else {
ss>>x;
if (ss.fail()) x=default_number; //use default if fail 
}
return x;
};

/*! \brief ask user for variable x with question in variable tx, present value of x is default
*/
template <class number>
number& Inquire(number &x, const char *tx){
return Inquire(x,tx,x);
};

/*! \brief boolean output to quesiton tx, there is also default option 
*/
bool Inquire(const char *tx, const bool default_bool) {
char xc;
cerr<<tx;
char default_char;
if (default_bool) default_char='y'; else default_char='n';
Inquire(xc,"(y/n)", default_char);
if (xc=='y') return true; else {if (xc=='n') return false; else return default_bool; }
} 

/*! 
\brief This allows to input data from a list of options in vector.
*/
template <typename data>
data Inquire(const std::vector<data> &offers) {
//for (int i=0;i<offers.size();i++) { 
//cerr<<offers[i]<<std::endl;
//}
cerr<<std::endl;
int i=0;
while (true) {
if (i==0) cerr<<"\nSelect your choice (type selection | <ret>=next | y=select) \n";
cerr<<offers[i]<<" ";
char c=cin.peek();
if (c!='\n') {
data entry;
cin>>entry;
for (int q=0;q<offers.size();q++) { 
if (entry==offers[q]) {
cerr<<"\nYour selection is: "<<entry<<std::endl; return entry;}
}
cerr<<"\nYour selection is: "<<offers[i]<<std::endl; return offers[i];
}
cin.ignore(); //this is return, so ignore it
i++;
i=(i%offers.size());
}
}

/*! 
\brief This allows to input data from a list of options in vector.
*/
template <typename containerx, typename data>
data InquireSelection(data &output, const containerx &offers) {
//for (int i=0;i<offers.size();i++) { 
//cerr<<offers[i]<<std::endl;
//}
cerr<<std::endl;
auto i=offers.begin();
while (true) {
if (i==offers.begin()) cerr<<"\nSelect your choice (type selection | <ret>=next | y=select) \n";
cerr<<(*i)<<" ";
char c=cin.peek();
if (c!='\n') {
data entry;
cin>>entry;
for (auto  q=offers.begin();q!=offers.end();++q) {
if (entry==(*q)) {
cerr<<"\nYour selection is: "<<entry<<std::endl; output=entry; return entry;}
}
cerr<<"\nYour selection is: "<<(*i)<<std::endl; output=(*i); return (*i);
}
cin.ignore(); //this is return, so ignore it
i++;
if (i==offers.end()) i=offers.begin();
}
}


}


#endif //__INQUIRE__CXX__
