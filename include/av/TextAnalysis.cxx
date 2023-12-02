/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/**
 \author Alexander Volya  <http://www.volya.net>
 \file TextAnalysis.cxx
 \brief Analysis of text configuration files
 */
/*
 02/23/07 fix fail bit at eof.
 DataStream.clear(); //clear fail bit
 06/27/09 use const *char
 01/19/10 CommentString using simple function of strspn and allow unlimited number of blank characters
 */
#ifndef TEXTANALYSIS__CXX__
#define TEXTANALYSIS__CXX__
#include<fstream>
#include <debug.h>
/*! \namespase av
 \brief non-templated library of functions and containers
 */
namespace av{
    using std::ifstream;
    using std::ios;
    using std::cout;
    ///\brief Characters identifying comments lines
    const char LocatePosition_Exclude_Char[] ="#!%/";
    /**\brief Check if string starts with a comment character, empty line is comment
     \return true-comment, false-not a comment
     */
    bool CommentString(const char *strFileLine ///<string
    ){
        int CharPos = std::strspn (strFileLine," "); //position of first character that is not space
        return (std::strspn (strFileLine+CharPos,LocatePosition_Exclude_Char)!=0); //if comment the the span of special characters is not zero
    };
    
    /** \brief Search for a first occurence of a character string in another string
     \return position in the line
     */
    unsigned int StringFind(const char* Sstring, ///<search string
                            const char* FullString ///< full string
    ) {
        const char terminate='\0';
        unsigned int pos=0;
        unsigned int scani=0;
        while (FullString[pos]!=terminate) {
            if (FullString[pos]==Sstring[0]) {
                do {
                    scani++;
                    if (Sstring[scani]==terminate) return pos; //success
                    if (FullString[pos+scani]==terminate)
                        return pos+scani; //no good
                }
                while (FullString[pos+scani]==Sstring[scani]);
                scani=0;
            }
            pos++;
        }
        return pos;//no good
    }
    
    /** \brief Skip Comments, new position is set after last comment line
     \return Number of comment lines, -1 is only coments
     */
    long SkipComments(std::ifstream &DataStream ///<file-stream
    ) {
        const int MaxLineLength=1024;
        char   strFileLine[MaxLineLength]; // One line buffer...
        unsigned long   ulLineCount=0;
        unsigned int   LinePos=0;
        long FilePosition;
        //std::ifstream DataStream;
        if (!DataStream.is_open()) FatalError("FileStringPositionFind:file is not open");
        if (DataStream){
            FilePosition=DataStream.tellg();
            while (DataStream.getline(strFileLine, MaxLineLength, '\n')){
                if (!(CommentString(strFileLine))) {
                    //      std::cerr<<strFileLine<<std::endl;
                    DataStream.seekg(FilePosition);
                    return ulLineCount;
                }
                ulLineCount++;
                FilePosition=DataStream.tellg();
            }
        }
        //DataStream.close();
        // DataStream.seekg(std::ios::end);
        DataStream.clear(); //clear fail bit
        
        return -1;
    }
    
    /** \brief Search in an open file stream a first occurence of a string.
     \return Position integer, if not found return -1.
     Assuming that the stream is open search first occurence of Sstring
     move position pointer to the location of the spring
     if string is not found return -1
     Comment lines are ignored.
     2/3/2016: Added a returnFlag argument that defaults to false. If set 
     to true, this will return the stream to wherever it was originally in 
     case the string was not found before returning -1.
     */
    long FileStringPositionFind(const char *Sstring, ///<search string
                                std::ifstream &DataStream, ///<file-stream
                                bool returnFlag = false///<Flag to control return in case of not found
    ) {
        const int MaxLineLength=1024;
        char   strFileLine[MaxLineLength]; // One line buffer...
        unsigned long   ulLineCount=0;
        unsigned int   LinePos=0;
        int LineLength;
        long FilePosition, startOfStream;
        //std::ifstream DataStream;
        if (!DataStream.is_open()) FatalError("FileStringPositionFind:file is not open");
        if (DataStream && returnFlag) startOfStream = DataStream.tellg();
        if (DataStream)
        {
            FilePosition=DataStream.tellg();
            while (DataStream.getline(strFileLine, MaxLineLength, '\n'))
            {

                /*get line discards last \n char, but inserts \0 as null terminating
                 the LineLength is counted with all chars */

                //ignore comment lines
                //********Check special chars************
                if (CommentString(strFileLine))
                {
                    FilePosition = DataStream.tellg();
                    continue;
                }
                //*********end special chars*************
                LinePos = StringFind(Sstring, strFileLine);
                if (strFileLine[LinePos] != '\0') {
                    DataStream.seekg(FilePosition + LinePos);
                    return FilePosition + LinePos;}
                //  else std::cout<<ulLineCount<<" "<<LinePos<<std::endl;
                
                //ignore empty line
                ulLineCount++;
                FilePosition=DataStream.tellg();
            }
        }
        //DataStream.close();
        // DataStream.seekg(std::ios::end);
        DataStream.clear(); //clear fail bit
        if (returnFlag) DataStream.seekg(startOfStream);
        return -1;
    }
    
    /**\brief Search for position with FileStringPositionFind() in file given with string name.
     */
    long FileStringPositionFind(const char * Sstring, const char * strFile) {
        std::ifstream DataStream;
        DataStream.open(strFile, std::ios::in|std::ios::binary);
        long FilePosition=FileStringPositionFind(Sstring, DataStream);
        DataStream.close();
        return FilePosition; //ulLineCount;
    }
    
    
    
    /**\brief Count lines in a text file, comments are ignored
     */
    unsigned long LineCount(const char * strFile) {
        const int MaxLineLength=1024;
        char   strFileLine[MaxLineLength]; // One line buffer...
        unsigned long   ulLineCount=0;
        int LineLength;
        std::ifstream DataStream;
        DataStream.open(strFile, std::ios::in|std::ios::binary);
        if (DataStream){
            while (DataStream.getline(strFileLine, MaxLineLength, '\n')){
                //LineLength=DataStream.gcount(); //count how many char were actually read
                /*get line discards last \n char, but inserts \0 as null terminating
                 the LineLength is counted with all chars */
                //std::cout<<strFileLine<<" "<<DataStream.gcount()<<std::endl;
                // cout<<int(strFileLine[0])<<'\n';
                // if (strFileLine[0]=='\0') continue; //empty line
                //ignore comment lines
                //if (strFileLine[0]=='!') continue;
                //if (strFileLine[0]=='#') continue;
                if (CommentString(strFileLine)) continue;
                //ignore empty line
                for (int i=0;strFileLine[i]!='\0';i++)
                    if (strFileLine[i]!=' ') {
                        ulLineCount++;
                        break; //preak from inner loop
                    }
            }
        }
        DataStream.close();
        return ulLineCount;
    }
    
    //note we add carrage return to protect from dos probelm 
#define IS_PUNC(c) ((c) == ' ' || (c) == ',' || (c)=='\r')
    
    /*** \brief Count number of words in the string. 
     */
    int WordCount(char *str)
    {
        int count = 0;
        char *s;
        
        // Skip leading punctuation and white-space
        while(IS_PUNC(*str))
            str++;
        s = str;
        
        while(*str)
        {
            while(*s && !IS_PUNC(*s))
                s++;
            if(str - s)
                count++;
            while(IS_PUNC(*s))
                s++;
            str = s;
        }
        
        return count;
    }
#undef IS_PUNC
    /**
     count word in current line of stream, starting from the current position
     */
    int WordCount(std::ifstream &DataStream) {
        if (!DataStream.is_open()) CriticalError("WordCount:file is not open");
        const int MaxLineLength=1024;
        char   strFileLine[MaxLineLength];
        if (!DataStream) CriticalError("WordCount:problem with filestream");
        long FilePosition=DataStream.tellg(); //save current position
        DataStream.getline(strFileLine, MaxLineLength, '\n');
        int words=WordCount(strFileLine);
        DataStream.clear(); //clear possible EOF
        DataStream.seekg(FilePosition); //return back to original position
        return words;
    }
} //end namespace
#endif

