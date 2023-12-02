/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file iofile.cxx
\brief String based checking file existance and reading from path
*/
#ifndef __IOFILE__CXX__
#define __IOFILE__CXX__
#include <debug.h>
using std::cerr;
using std::endl;
#include <fstream> //this must be included
using std::fstream;
#include <tav/strings.cxx>
#include <vector>
//#include <string>
//using std::string;
namespace av{
/// Check File existance
bool CheckFileExistance(const string &filename)
{
  std::ifstream ifile(filename.c_str());
  bool notfail=(!ifile.fail()); //check status 
  if (notfail) {
  ifile.get(); //get character for certainty
  notfail=ifile.good(); //confirm that everything is fine
  ifile.close(); //close if open was successful
  }
  return notfail;
}

/*
//it is possible to use the following funciton 
void TokenizeString(const string& str,
                      stl::vector<strings>& tokens,
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
*/
/** \brief Open file by file name and try string of path
 */
string FullFileName(
  const string &str, ///< [i/0] input: list of paths to try, output: succesful name MUST BE LARGE ENOUG!!
  const string &filename  ///< [i] short file name
)
{
  if (CheckFileExistance(filename)) return filename;
 // if (std::strlen(str)==0) CriticalError("FullFileName: File cannot be found");
  string separators=";:";
  std::vector<string> paths;
  tav::TokenizeString(str,paths, separators); //break up string
  for (int i=0;i<paths.size();i++) {
  string fullfilename=paths[i];
  if (paths[i].length()>0) if (paths[i][paths[i].length()-1]!='/') fullfilename+="/"; //add "/" at end if needed
  fullfilename+=filename;
  if (CheckFileExistance(fullfilename)) return fullfilename;
  }
  
   //we have an error
  CriticalError("FullFileName: File "<<filename<<" cannot be found \n within path ");
  return str;
}

/** \brief Code returns a file name based on short name and pathvariable
returns c-string to be compatable with file.open(name)
*/
string FileName(
const string &shortname, ///< [i] Short name of the file
const string &envname, ///< [i] Path variable
const string &extension="" /// [i] Possible extension
) {
string pathstring;
  if (std::getenv(envname.c_str())!=NULL) pathstring=std::getenv(envname.c_str());
  tav::RemoveNonPrint(pathstring); //remove non-printable characters
  try {
  return av::FullFileName(pathstring, shortname);
  }
  catch (...) {
  //prepare name with extension
  string longertmp=shortname;
  if (extension[0]!='.') longertmp+="."; //add "." if needed
  longertmp+=extension;
  //read again path string
  // if (std::getenv(envname)!=NULL) std::strcpy(pathstring, std::getenv(envname));
  try {return av::FullFileName(pathstring, longertmp);}
  catch (...) {
  cerr<<"File not found: "<<shortname<<", check path: "<<envname<<"\n";
//  cerr<<"---------------------File Not Found --------------------------"<<endl;
//  cerr<<"test name: "<<shortname<<"\n";
//  cerr<<"test name: "<<longertmp<<"\n";
//  cerr<<"path environment variable: "<<envname<<"\n";
// // cerr<<"search path: \n"<<pathstring<<endl;
//  cerr<<"directories searched\n";
//  string separators=";:";
//  std::vector<string> paths;
//  tav::TokenizeString(pathstring,paths, separators); //break up string
//  for (int i=0;i<paths.size();i++) cerr<<paths[i]<<endl;
	//print ouput error
	//if (std::getenv(envname)!=NULL) std::strcpy(pathstring, std::getenv(envname));
	CriticalError("FullFile: File cannot be found ");
    }
  }
}


/** \brief Relative path to a file retrun with trailing slash:  "/"
 */
string FileSearchPath(
  const string &str, ///< [i] path variable
  const string &filename  ///< [i] short file name
)
{

  if (CheckFileExistance(filename)) return "./"; //if file is here return relative dir
 // if (std::strlen(str)==0) CriticalError("FullFileName: File cannot be found");
 string retpath; //return path
  string separators=";:";
  std::vector<string> paths;
  tav::TokenizeString(str,paths, separators); //break up string
  for (int i=0;i<paths.size();i++) {
  string fullfilename=paths[i];
  if (paths[i].length()>0) if (paths[i][paths[i].length()-1]!='/') fullfilename+="/"; //add "/" at end if needed
  //fullfilename+=filename;
  if (CheckFileExistance(fullfilename+filename)) return fullfilename;
  }
  
   //we have an error
  CriticalError("FilePath: file "<<filename<<" cannot be found \n within path ");
  return str;
}

string FilePath(
  const string &envname, ///<  path variable name
  const string &filename  ///< [i] short file name
)
{
string pathstring;
if (std::getenv(envname.c_str())!=NULL) pathstring=std::getenv(envname.c_str());
return FileSearchPath(pathstring,filename);
}
} //end namespace av
#endif //__IOFILE__CXX__

