/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file strings.cxx
\ingroup gp_av
*/
#ifndef __STRINGS__CXX__
#define __STRINGS__CXX__
#include <cctype> //or use <locale>
using std::isprint;
#include <string>
using std::string;
#include <algorithm>
using std::remove_if;
namespace tav{
 bool isnotprint(int c) {return (! std::isprint(c));}
//Remove all nonprint characters
int RemoveNonPrint(string &s) { s.erase(std::remove_if(s.begin(), s.end(), (isnotprint)), s.end()); return 0;}  
//Clean dublicate spaces
int CleanSpaces(string &mystring) {
for (int c=0;c<mystring.length();)  if (mystring[c]==' ')  mystring.erase(c,1);  else break;   //note no need to increment
for (int c=mystring.length()-1 ;c>=0 ;c--)  if (mystring[c]==' ')  mystring.erase(c,1);  else break;
return 0;
}
/**! \breif Split string into tokens(works) and save them in svector, any class that allows pushback operator
*/
template <typename svector> //std::vector<string>
void TokenizeString(const string& str,
                      svector& tokens,
                      const string& delimiters = " ")
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}
/**! \brief obtain ith token(word) from a string, empty string if not found 
*/
string StringToken(int i, const string& str,
                      const string& delimiters = " ")
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);
    int ii=0; 
    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        //tokens.push_back(str.substr(lastPos, pos - lastPos));
		if (i==ii) return str.substr(lastPos, pos - lastPos);
		ii++;
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
	return string(); //return empty if not found
}
} //end namespace tav
#endif //__STRINGS__CXX__

