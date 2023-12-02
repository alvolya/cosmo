/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*!
 \author Alexander Volya  <http://www.volya.net>
 \file stensor.cxx
 \brief STensor class and associated functions
 \example stensor_test.cpp
 \brief An example of basic STensor usage, also Read/Write
 */
/** \example stensor_access_test.cpp
 This is an example
 */
/*
 07/16/09 operator >> (for read tensor is introduced) it is basic and requires dimension first
 07/16/08 Read and Write non-member introduced
 07/15/08 cnumber(0.0) is changed to cnumber(0) in make_pair
 07/07/2007
 New Version of STensor, rank can be a variable, no elem, use map as base class
 07/10/2007 new variable mem_size is included, copy memory array.
 memcpy can outperform significantly on large arrays, however comparison is very slow. strcmp!!
 10/22/2007 Operator *=
 10/23/2007 Operator +=
 12/29/2010 Operator *(double)
 12/29/2010 Operator /(double)
 12/29/2010 Operator -, fix with constant
 12/29/2010 Operator +
 09/17/2016 Read write modified, allows to read rank, add SetRank function
 */
#ifndef __STENSOR__CXX__
#define __STENSOR__CXX__
#include <debug.h>
//debug calls iostream
//use for FatalError
#include <map>
#include <cstdarg>
//#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>    // std::swap
//#include <limits> // std::numeric_limits<double>::epsilon()
#include <tav/basic.cxx>
using tav::Sqr;
using tav::Abs;
#include <rw/basic-rw.cxx> //need TestFile
using rw::TestFile;
using namespace std;
namespace tav {
    //    typedef unsigned char base_integer; ///< integer for denoting rank [unsigned char]
    typedef uint16_t base_integer;
    /*! \brief Helper comparison class for STensor
     
     Provides a comparison (less) function operator()
     */
    template <typename ibase>
    class ArrayComparison
    {
    public:
        base_integer rank_C; ///< rank or the array
        
        /// Default constructor assumes rank=1
        ArrayComparison(base_integer _rank=1) : rank_C(_rank) {}
        
        /*!
         \brief compare two arrays given by pointers, true if *s1 < *s2
         */
        bool operator()(const ibase* s1, const ibase* s2) const
        { //cerr<<" rank here "<<rank<<endl;
            for (base_integer i=0;i<rank_C;i++)
                if (s1[i]!=s2[i]) return (s1[i]<s2[i]);
            return false;
        }
        void SetRank(base_integer _rank){rank_C = _rank;}
        
    };
    
    
    
    
    /*! \class STensor
     \brief STensor is a multidimensional sparse matrix/tensor
     
     STensor is a class that uses std::map as its base. It maps an array of arguments [domain] onto a output variable [range]. The purpose of this class is to provide
     - User defined comparison #ArrayComparison
     - Dynamic rank, which can be defined/change during runtime
     - Creation and deletion of new elements with new/delete of array arguments
     - Sorting, search, iterators are inherited from std::map
     */
    template <typename ibase, typename cnumber>
    class STensor : public std::map<ibase*, cnumber, ArrayComparison<ibase> > {
    private:
        
    public:
        base_integer rank; //we can not override this
        
        int elem_size; ///< number of bytes in an element
        typedef std::map<ibase*, cnumber, ArrayComparison<ibase> > Map;
        typedef typename Map::iterator iterator;
        typedef typename Map::size_type size_type; //size_type can be statically defined for saving and reading
        //try this as default constructor
        STensor() : Map(0), rank(0),
				    elem_size(0) {};
        
        //constructor inheritance
        STensor(base_integer _rank) : Map(_rank), rank(_rank),
				    elem_size(_rank*sizeof(ibase)) {};
        //copy constructor
        STensor(const STensor &T) : Map(T.rank), rank(T.rank),
				    elem_size(T.rank*sizeof(ibase))
        {
//            cerr<<"Copy constructor called "<<endl;
            ibase *dummy;
            typename STensor<ibase,cnumber>::const_iterator iter;
            for(iter = T.begin(); iter != T.end(); iter++) {
                dummy=new ibase[rank]; //create new empty array
                memcpy(dummy,iter->first,elem_size); //copy this array
                //insert at the end.
                Map::insert((*this).end(),make_pair(dummy,iter->second));
            }
        };
        //destructor
        ~STensor(void) {
            //      std::cout << "Entering destructor" << std::endl;
            for(iterator iter = Map::begin(); iter != Map::end(); iter++)
                delete [] (*iter).first;
        };
        ///\brief function to clear STensor, similar to std::map::clear()
        void clear(void) {
            for(iterator iter = Map::begin(); iter != Map::end(); iter++)
                delete [] (*iter).first;
            Map::clear();
        };
        ///\brief function to clear STensor, similar to std::map::erase(), fix to return erase
        void erase(iterator iter) { delete [] (*iter).first; Map::erase(iter); };
        
   //change rank of tensor, if rank is the same then tensor remains unchanged else cleared
        base_integer SetRank(base_integer _rank) {
        if (_rank!=rank) {
        this->clear();
        tav::STensor<ibase,cnumber> tq(_rank);
            *this = tq;
        }
        return rank;
        }
        
        /*****************INSERT AND MORE********************/
        /*! \brief element access operator.
         Equivalent to map::operator [],
         but assures creation, deep copy of the argument array and inception in map.
         Argument is constant, the return element is unchanged
         \note Often it is necessary to force copy so that the original pointer is not lost then c++ cast must be done.
         const_cast<const ibase*>(a)
         */
        cnumber& operator [] (const ibase* a) {
            
            /*
             no second locate but we still need copy before
             T& operator[] ( const key_type& x ); is equal to
             (*((this->insert(make_pair(x,T()))).first)).second
             inx is a result of incert; if new lelement
             
             */
            ibase *dummy=new ibase[rank]; //create new empty array
            //for (base_integer i=0;i<rank;i++) dummy[i]=a[i];
            memcpy(dummy,a,elem_size); //copy this array
            pair<iterator,bool> inx = Map::insert(make_pair(dummy,cnumber(0)));
            if (!inx.second) {//not new
                delete [] dummy;
            }
            return (*(inx.first)).second;
            
        };
        
        pair<iterator,bool> insert (const ibase* a) {
            
            /*
             no second locate but we still need copy before
             T& operator[] ( const key_type& x ); is equal to
             (*((this->insert(make_pair(x,T()))).first)).second
             inx is a result of incert; if new lelement
             
             */
            ibase *dummy=new ibase[rank]; //create new empty array
            //for (base_integer i=0;i<rank;i++) dummy[i]=a[i];
            memcpy(dummy,a,elem_size); //copy this array
            pair<iterator,bool> inx = Map::insert(make_pair(dummy,cnumber()));
            if (!inx.second) {//not new
                delete [] dummy;
            }
            return inx;
            
        };
        
        
        pair<iterator,bool> insert (ibase*& a) {
            pair<iterator,bool> inx = Map::insert(make_pair(a,cnumber()));
            if (inx.second) a=new ibase[rank]; //create new empty array
            return inx;
        };
        
        
        /*! \brief operator [], not constant, reference to pointer.
         Oparator is faster since it only inserts the string and creates new space
         when incertion is successful.
         Returned pointer a is EMPTY, initial data is not preserved !
         */
        cnumber& operator [] (ibase*& a) {
            pair<iterator,bool> inx = Map::insert(make_pair(a,cnumber(0)));
            if (inx.second) a=new ibase[rank]; //create new empty array
            return (*(inx.first)).second;
        };
        
        //this only is needed because of other find below
        //  inline iterator find(ibase* a) {return Map::find(a);};
        //my find goes with the capital
        ///\brief cstdarg from of access, returns iterator to map.
        iterator Find(ibase n1,...) {
            va_list ap;
            va_start(ap,n1); //initialize list
            ibase *dmn=new ibase [rank]; //create an array for data
            dmn[0]=n1;
            for (int i=1;i<rank;i++) dmn[i]=ibase(va_arg(ap,int));
            //note variable argument requires integer
            va_end(ap);
            iterator iit=find(dmn);
            delete [] dmn;
            return iit;
        }
        
        
        //this can be simplified because of operator [] on map acts as needed
        //!. \brief cstdarg access to an element, returns output type
        cnumber & operator () (ibase n1, ...) {
            va_list ap;
            va_start(ap,n1); //initialize list
            ibase *dmn=new ibase [rank]; //create an array for data
            dmn[0]=n1;
            for (int i=1;i<rank;i++) dmn[i]=ibase(va_arg(ap,int));
            //note variable argument requires integer
            va_end(ap);
            iterator iit=Map::find(dmn);
            if (iit==Map::end()) {
                //did not find, assignment only possible we keep the array dmn
                return Map::operator [] (dmn);
            }
            else {
                delete [] dmn; //in this case just return the value and delete dmn
                return (*iit).second;
            }
        }
        /*! \brief multiply STensor by a number
         \return reference to the resulting STensor;
         */
        template <typename qnumber>
        STensor&  operator *= (const qnumber &T) {
            for (iterator iter=Map::begin(); iter !=Map::end(); iter++) ((*iter).second)*=T;
            return *this;
        };
        
        /*! \brief multiply STensor by a number
         \return constant  STensor, note there is copy done twice;
         */
        template <typename qnumber>
        STensor  operator* (const qnumber &T) const
        {
            return STensor<ibase,cnumber>(*this)*=T;
        };
        
        
        
        
        /*! \brief divide STensor by a number
         \return reference to the resulting STensor;
         */
        template <typename qnumber>
        STensor&  operator /= (const qnumber &T) {
            for (iterator iter=Map::begin(); iter !=Map::end(); iter++) ((*iter).second)/=T;
            return *this;
        };
        /*! \brief divide STensor by a number
         \return constant  STensor, note there is copy done twice;
         */
        template <typename qnumber>
        STensor  operator/ (const qnumber &T) const
        {
            return STensor<ibase,cnumber>(*this)/=T;
        };
        
        /*! \brief Add STensors
         \return reference to the resulting STensor;
         */
        template <typename qnumber>
        STensor&  operator += (const STensor<ibase,qnumber> &T) {
            if (rank!=T.rank) FatalError("operator +=STensor<>: ranks are not equal");
            typename STensor<ibase,qnumber>::const_iterator iter;
            //typename STensor<ibase,cnumber>::iterator iter;
            for(iter = T.begin(); iter != T.end(); iter++)
                (*this)[(iter->first)]+=cnumber(iter->second);
            ///\todo speed can be improved, use position
            return *this;
        };
        
        /*! \brief Add STensors
         \return reference to the resulting STensor;
         
         */
        template <typename qnumber>
        STensor  operator + (const STensor<ibase,qnumber> &T) const {
            return STensor<ibase,cnumber>(*this)+=T;	;
        };
        
        /*! \brief Subtracting STensors
         \return reference to the resulting STensor;
         \note can not use "const STensor &T"
         */
        template <typename qnumber>
        STensor&  operator -= (const STensor<ibase,qnumber> &T) {
            if (rank!=T.rank) FatalError("operator +=STensor<>: ranks are not equal");
            typename STensor<ibase,qnumber>::const_iterator iter;
            //typename STensor<ibase,cnumber>::iterator iter;
            for(iter = T.begin(); iter != T.end(); iter++)
                (*this)[(iter->first)]-=cnumber(iter->second);
            ///\todo speed can be improved, use position
            return *this;
        };
        
        /*! \brief subtracting STensors
         \return reference to the resulting STensor;
         \note can not use "const STensor &T"
         */
        template <typename qnumber>
        STensor  operator - (const STensor<ibase,qnumber> &T) const {
            return STensor<ibase,cnumber>(*this)-=T;	;
        };
        
        
        /*! \brief copy STensors
         \return reference to the resulting STensor;
         \note NOT well tested !!!
         */
        STensor&  operator = (const STensor &T)
        {
//            cerr<<" equal called "<<endl;
            this->clear(); //remove everything
            this->rank=T.rank;
            this->elem_size=T.rank*sizeof(ibase);
            //Map tmp(rank);
            //Map(*this)=tmp;
            
            Map temp(T.rank);
            
            for (auto it = T.begin(); it!= T.end(); ++it)
            {
                ibase *dummy=new ibase[T.rank]; //create new empty array
                memcpy(dummy,it->first,elem_size); //copy this array
                temp[dummy] = it->second;
            }
            
            Map* mapptr = &temp;
            Map* tempthis = this;
            std::swap(*tempthis,*mapptr);
            
            ///\todo speed can be improved, use position
            
            return *this;
        };  /**********************************************/
        
        
        int Write(const char *filename);
        int Write(std::ofstream&);
        int Read(const char *filename);
        int Read(std::ifstream&);
        
    };
    
    
    //! Save to stream (binary)
    template <typename ibase, typename cnumber>
    int STensor<ibase, cnumber>::Write(std::ofstream &outf) {
        
        outf.write((char*)(&rank),sizeof(base_integer)); //save rank first
        size_type wsize=Map::size();
        outf.write((char*)(&(wsize)),sizeof(size_type)); //save size second
        int sz=rank*sizeof(ibase);
        for (iterator iter=Map::begin(); iter !=Map::end(); iter++) {
            outf.write((char*)((iter->first)),sz);
            outf.write((char*)(&(iter->second)),sizeof(cnumber));
        }
        return 0;
    };
    
    //! Write to stream (binary), apply write function to all elements
    template <typename ibase, typename cnumber>
    int Write(STensor<ibase, cnumber> &T, std::ofstream &outf) {
        outf.write((char*)(&(T.rank)),sizeof(base_integer)); //save rank first
        auto wsize=T.size();
        outf.write((char*)(&(wsize)),sizeof(typename STensor<ibase, cnumber>::size_type));
        int sz=T.rank*sizeof(ibase);
        typename STensor<ibase, cnumber>::iterator iter;
        for (iter=T.begin(); iter !=T.end(); iter++) {
            //for (base_integer i=0;i<T.rank;i++) tav::Write(iter->first[i],outf);
            outf.write((char*)((iter->first)),sz); //use above if structure is more complex
            tav::Write(iter->second, outf);
            //outf.write((char*)(&(iter->second)),sizeof(cnumber));
        }
        return 0;
    };
    
    //! Save to file (binary)
    template <typename ibase, typename cnumber>
    int Write(STensor<ibase, cnumber> &T, const char* filename) {
        std::ofstream outf;
        outf.open(filename, std::ios::binary);
        int ret=Write(T,outf);
        outf.close();
        return ret;
    };
    
    //! Save to file (binary)
    template <typename ibase, typename cnumber>
    int STensor<ibase, cnumber>::Write(const char* filename) {
        std::ofstream outf;
        outf.open(filename, std::ios::binary);
        int ret=Write(outf);
        outf.close();
        return ret;
    };
    
    //! Read from stream (binary)
    template <typename ibase, typename cnumber>
    int STensor<ibase, cnumber>::Read(std::ifstream &inf) {
        this->clear(); //remove everything
        base_integer inrank; //read value of rank
        inf.read((char*)(&inrank), sizeof(base_integer)); //read rank
        //cerr<<"Read rank "<<inrank<<endl;
        size_type insize, readsize=0; //read size
        inf.read((char*)(&insize), sizeof(size_type)); //read size
        //cerr<<"Saved size "<<insize<<endl;
        //if new rank is not the same make changes
        
//        if (this->rank != inrank) {
//            tav::STensor<ibase,cnumber> tq(inrank);
//            *this = tq;//equality fills this with empty tensor of correct rank.
//        }
        this->SetRank(inrank);
        
        //proceed reading
        int sz=rank*sizeof(ibase);
        ibase *dummy_first;
        cnumber dummy_second;
        char _ctmp=inf.peek();
        //char _ctmp;
        while (!inf.eof())
        {
            dummy_first=new ibase[rank]; //each run make new array
            inf.read((char*)(dummy_first),sz); //we read first element in dummy array;
            inf.read((char*)(&dummy_second), sizeof(cnumber));
            Map::operator [] (dummy_first)=dummy_second; //insert here
            readsize++;
            //cerr<<(*this).size()<<" "<<dummy_second<<" "<<dummy_first[0]<<" "<<dummy_first[1]<<" "<<dummy_first[2]<<" "<<dummy_first[3]<<endl;
            _ctmp=inf.peek();
        }
        //cerr<<"Read size "<<readsize<<endl;
        if (readsize!=insize)  FatalError("Stensor.Read(stream) size mismatch");
        return 0;
    };
    
    
    //! Read from stream (binary) with read Format
    template <typename ibase, typename cnumber>
    int Read(STensor<ibase, cnumber> &T, std::ifstream &inf) {
        T.clear();
        base_integer inrank; //read value of rank
        inf.read((char*)(&inrank), sizeof(base_integer)); //read rank
        typename STensor<ibase, cnumber>::size_type insize, readsize=0; //read size
        inf.read((char*)(&insize), sizeof(STensor<ibase, cnumber>::size_type)); //read size
        
        //if new rank is not the same make changes
//        if (T.rank!=inrank) {
//            T=STensor<ibase, cnumber> (inrank); //create new variable
//        }
        T.SetRank(inrank);
        
        int sz=T.rank*sizeof(ibase);
        ibase *dummy_first;
        cnumber dummy_second;
        char _ctmp=inf.peek();
        //char _ctmp;
        while (!inf.eof() )
        {
            dummy_first=new ibase[T.rank]; //each run make new array
            //for (base_integer i=0;i<T.rank;i++) tav::Read(dummy_first[i],inf);
            inf.read((char*)(dummy_first),sz); //we read first element in dummy array;
            // inf.read((char*)(&dummy_second), sizeof(cnumber));
            rw::Read(T[dummy_first],inf);
            readsize++;
            //Map::operator [] (dummy_first)=dummy_second; //incert here
            //T[dummy_first]=dummy_second;
            _ctmp=inf.peek();
        }
        if (readsize!=insize)  FatalError("Read(STensor, stream) size mismatch");
        return 0;
    };
    
    //! Read from file (binary)
    template <typename ibase, typename cnumber>
    int STensor<ibase, cnumber>::Read(const char* filename) {
        if (!TestFile(filename))
            FatalError("Stensor.Read(file) need file to exist");
        std::ifstream inf;
        inf.open(filename, std::ios::binary);
        int ret=Read(inf);
        inf.close();
        return ret;
    };
    
    //! Read from file (binary)
    template <typename ibase, typename cnumber>
    int Read(STensor<ibase, cnumber> &T, const char* filename) {
        if (!TestFile(filename))
            FatalError("Read(STensor,file) need file to exist");
        std::ifstream inf;
        inf.open(filename, std::ios::binary);
        int ret=Read(T,inf);
        inf.close();
        return ret;
    };
    
    template<typename ibase, typename cnumber>
    std::ostream& operator << (std::ostream &stream,
                               STensor<ibase, cnumber> &T)
    {
        typename STensor<ibase, cnumber>::iterator iter;
        stream << T.size()<<std::endl;
        for(iter = T.begin(); iter != T.end(); iter++)
        { for (base_integer i=0;i<T.rank;i++) stream<<int((*iter).first[i])<<" ";
            stream<<(*iter).second<<std::endl;
        }
        return stream;
    }
    
    template<typename ibase, typename cnumber>
    std::istream& operator >> (std::istream &DataStream,
                               STensor<ibase, cnumber> &T)
    {
        T.clear();
        int size;
        DataStream>>size; //read array size
        ibase *a=new ibase [T.rank];
        cnumber x; //variables to read
        for (int i=0;i<size;i++) {
            for (int i=0;i<T.rank;i++) DataStream>>a[i];
            DataStream>>x;
            if (!DataStream.good()) {
                if (DataStream.fail()) FatalError("stream>>Stensor: Loss of stream integrity while reading \n");
                if (DataStream.eof()) FatalError("stream>>STensor: End-of-file reached while reading \n");
            }
            //manipulate
            T[a]=x; //incert element
        }
        
        delete [] a;
        return DataStream;
    }
    
    /*
     //The following read should be considered with sstring (stream strings)
     #include <av/TextAnalysis.cxx> //use WordCount and CommentString
     using av::WordCount;
     using av::CommentString;
     template<typename ibase, typename cnumber>
     std::istream& operator >> (std::istream &DataStream,
     STensor<ibase, cnumber> &T)
     {
     T.clear();
     const int MaxLineLength=1024;
     char   strFileLine[MaxLineLength]; // One line buffer...
   	 ibase *a=new ibase [T.rank];
     cnumber x; //variables to read
     int words;
     if (DataStream){
     long FilePosition=DataStream.tellg();
     while (DataStream.getline(strFileLine, MaxLineLength, '\n')){
     cerr<<"We got string "<<strFileLine<<endl;
     if (CommentString(strFileLine)) 
     { FilePosition=DataStream.tellg(); continue;}
     //-------------we do stuff-----------------
     words=WordCount(strFileLine);
     cerr<<"number of workds is "<<words<<endl;
     if (words!=T.rank+1) 
     { cerr<<"stop read"<<endl; 
     FilePosition=DataStream.tellg(); break;}
     //cout<<"line ="<<ulLineCount<<" words="<<words<<endl;
     DataStream.clear(); //clear possible EOF
     DataStream.seekg(FilePosition);
     //read and store data
     for (int i=0;i<words-1;i++) DataStream>>a[i];
     DataStream>>x; 
     //manipulate
     T[a]=x; //incert element
     DataStream.getline(strFileLine, MaxLineLength, '\n');
     //----------------end stuff
     FilePosition=DataStream.tellg();
     }
     }
     delete [] a;
     return DataStream;
     }
     
     */
    
    /*! \fn Project
     \brief View tensors as vectors and project one out of the second one
     \return projection coefficient
     
     \f[
     {\vec a}\rightarrow {\vec a}-c {\vec x}\,,\quad c=\frac{{\vec a}\cdot {\vec x}}{{\vec x}\cdot {\vec x}}
     \f]
     */
    template<typename ibase, typename cnumber>
    cnumber Project(
                    tav::STensor<ibase, cnumber> &VSP, ///<[input/output] tensor to be modified  
                    tav::STensor<ibase, cnumber> &PPP ///<[input] projected out component of interaction 
    ) {  
        cnumber ppnorm=0.0;
        cnumber pproduct=0.0;
        if (PPP.empty()) FatalError("Project: the STensor is empty");
        typename tav::STensor<ibase, cnumber>::iterator ip;
        //compute norm of projecting operator and product
        for (ip=PPP.begin();ip!=PPP.end();ip++) {
            ppnorm+=tav::Sqr((*ip).second);
            pproduct+=(*ip).second*VSP[ip->first]; //if empty this will create new zero-initilized entry
        }
        if (tav::Abs(ppnorm)<1E-10) FatalError("Project: the STensor is zero");
        cnumber component=pproduct/ppnorm;	
        for (ip=PPP.begin();ip!=PPP.end();ip++) {
            VSP[ip->first]-=component*(ip->second);
        }	
        
        return component;
    }
    
    /*! \fn Sqr
     \brief Square tensor
     */
    template<typename ibase, typename cnumber>
    cnumber Sqr(
                const tav::STensor<ibase, cnumber> &PPP ///<[input] tensor to be modified
    ) {  
        cnumber ppnorm=0.0;
        //typename tav::STensor<ibase, cnumber>::iterator ip;
        //compute norm of projecting operator and product
        for (auto ip=PPP.begin();ip!=PPP.end();ip++)
        ppnorm+=(ip->second)*(ip->second);
        return ppnorm;
    }
    
    
    /*! \brief Strips operator of Zero elements */
    template<typename ibase, typename cnumber>
    int StripZeros1(tav::STensor<ibase, cnumber> &input,double tol=1E-10)
    {
        tav::STensor<ibase, cnumber> temp(input.rank);
        std::swap(input,temp);
        input.clear(); //just to make sure
        for (auto it = temp.begin(); it != temp.end(); it++)
            if (std::abs(it->second) > tol) input[it->first]=it->second;
        return input.size();
    }
    
    template<typename ibase, typename cnumber>
    int StripZeros(tav::STensor<ibase, cnumber> &input,double tol=1E-10)
    {
        for (auto it = input.begin(); it != input.end();)
            if (std::abs(it->second) <= tol) {
                input.erase(it++); }
            else
            {
                ++it;
            }
        return input.size();
    }
    
    
    
} //end namespace tav



#endif

