/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file tensor.cxx
\brief Multidimensional tensor class. 

The class defines a storage container of multidimensional array (called tensor)
the dimensionality is given by the rank. 
- tav::Tensor<1,double> A(10); defines a simple array, A can be used as pointer. This is equivalent double *A=new double[10]
- tav::Tensor<2,double> B(5,4); defines two dimensional array B, equivalent
double **B= new double* [5]; for (int i=0;i<5;i++) B[i]=new double [4]; 
Access is usual B[i][j]. Destruction is controlled automatically. 
- Not generally OpenMP safe, but some steps have been taken. 

\example Tensor
Tensor class usage example 
\code
#include <iostream>
#include <fstream>
//#include <complex>
#include "tensor.cxx"
using namespace std;
int main(void) {
  tav::Tensor<1, double> CC(4), DD(5); //define CC[4] and DD[5] one dim arrays
    CC=7.0;  // set all values in array to be equal 7.0 
     cout<<CC[1]<<" "<<DD<<endl; //print CC[1] and entire array DD
     Swap(CC,DD); //swap arrays
      cout<<CC[1]<<" "<<DD[1]<<endl;
   tav::Tensor<3,double> Y(2,3,4); //define Y[2][3][4] 
   Y[1][0][2]=2.7; //assigns the element to 2.7
   return 1;
}
\endcode
 
 
 
 
 

 */
/*
04/12/2005 this is ansi C++ complient version of tensor global_pointer 
is introduced to pass constructor invormation 
operator equal for unequal dimensions will fail
09/23/2005 include Tensor<rank..> op number 
where op is += -= = *= and so on. 
    Resize function is available for 1D tensor
    1D read will resize appropriately
    ND resize is functioning
    ND read will resize appropriately
    constructor from file is written
    introduced global_zero array for default construction with zeros 
    (max rank 16)
    Extra function Full Resize is introduced 
    constactor Tensor(Tensor) is built
    operator = is fixed
09/30/2005 initilization of global_zero fixed pointer must go to the end
10/03/2005 there seems to be a problem with creation from file    
10/14/2005 include Read(tensor, file) and Write operators these are 
    by-element read-write operatons, should be compativle with .Read/Write
    fixed close stream error in Read/Write
10/17/2005 Operation Swap =Swap(dim)+ Swap(elem); works for arbitrary rank
10/18/2005 Expand function is included, similar to Resize but keep elements
10/19/2005 cast operator Tensor<1,cnumber> to cnumber*
10/20/2005 The problem of creating tensor from file is fixed, 
    elem=NULL required
12/17/2005 FullExpand INEFFICIENT!
    rank=1 Full expand from pointer
    FullExpand for any rank with pointer
    FullExpand for vararg
12/18/2005 Member function Dim(*int) return array of dimensions(rectangle type)
03/09/2006 problem with = operator this version works but noes not work on radon
03/16/2006 operator >> (istream, Tensor<1(rank),rank>) is introduced  will not change dimensions
04/18/2006 class changed to typename, forward declaration is introduced
04/19/2006 operator= fixed again, no need for Tensor::operator= in prototypes
04/06/2007 intruduce TimesVector operator for left and right
           connjugate multiplication??? | 
07/16/2007 all Read/Write now uses const char* 
09/26/2007 operator [] is removed, we use casting to a pointer 
elem can be private?
05/11/2008 error handling included,
constractor Tensor() is needed only for new[] , introduced initialize finction
07/16/2008 In copy constructor Tensor(const Tensor&) need const
 */
#ifndef __TENSOR_CXX__
#define __TENSOR_CXX__
//#include <iostream>
#include <fstream>
#include <cstdarg>
#include <debug.h>
//debug calls iostream
#include <tav/basic.cxx>  
/*dependence: Read/Write, Min, Swap */
using std::cerr;
using std::endl;
typedef  unsigned long long i_number ;///< Integer variable to describe size of each dimension
namespace tav{


/// global array of zeros, used in initializations
  static i_number global_zero[32]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; 
/// global pointer to array (end) of zeros, used in initializations
  i_number *global_pointer=(global_zero+31); //this golobal pointer helps control dimensions
    // Forward declarations
    //template<int rank, typename cnumber> class Tensor;
  /*! \class Tensor
    \brief Multidimensional data structure arrasy 
  */
    template <int rank, typename cnumber>
  class Tensor
  {
  private:
  /*! \brief initializing function
  */
       void initialize(i_number* dmn){
	    #pragma omp critical 
    {
    global_pointer=&dmn[rank-1];
    //this is first call and we store global pointer to last dimension
    dim=(global_pointer-rank+1)[0];
    //elem=new Tensor<rank-1,cnumber>[dim](dmn+1);
    try {elem=new Tensor<rank-1,cnumber>[dim];}
//process exception possibility
   catch (std::exception& e)
    {
	cerr<<"Exception: constructor Tensor(";
	for (int ii=0;ii<rank-1;ii++) cerr<<dmn[ii]<<",";
	cerr<<dmn[rank-1]<<")"<<endl;
    FatalError(e.what());
    }
	  global_pointer=(global_zero+31);
	 }//end omp critical region
	   }; //min constructor
	   
	   
  public: 
      i_number dim; ///< dimension top structure 
      Tensor<rank-1,cnumber>  *elem; ///< pointer to lower rank substructure
    //Constructors  
    //! default constructor set all dim=0 
    Tensor(); 
    //! constructor from array  dim[] 
    Tensor(i_number*); //min constructor
    //! constructor from list of arguments (dim1,dim2,...)  
    Tensor(i_number d1, ...); //cstdarg function multiple artuments
    //! construct tensor from file    
    Tensor(const char *filename) {elem=NULL; Read(filename);}; //reading from file
   //! construct tensor from another tensor    
    Tensor(const Tensor<rank,cnumber>& ); 
    ~Tensor(void); ///<destructor
   //! cast tensor to pointer for a lower rank structure, which allows [] access
    inline operator Tensor<rank-1,cnumber>* (void) const {return elem;};
     inline i_number size(void) {return dim;}  ///< size function for compatibility
  //! deep copy  identical temsor  
      Tensor<rank,cnumber>& operator=(const Tensor<rank, cnumber>& T); 
  //! deep copy casting tensor
    template<typename qnumber>
    Tensor<rank,cnumber>&  operator=(const Tensor<rank, qnumber>& T); 
  


 //Operators Tensor A+=number 
 template <typename qnumber>
    Tensor<rank,cnumber>&  operator *= (const qnumber T); 
 template <typename qnumber>
    Tensor<rank,cnumber>&  operator /= (const qnumber T);
 template <typename qnumber> 
    Tensor<rank,cnumber>&  operator = (const qnumber T);
 template <typename qnumber>
    Tensor<rank,cnumber>&  operator += (const qnumber T);
 template <typename qnumber>
    Tensor<rank,cnumber>&  operator -= (const qnumber T);



//! Write tensor to file (binary)
    int Write(const char *filename);
//! Write tensor to stream (binary)
    int Write(std::ofstream&);
//! Read tensor from file (binary)
    int Read(const char *filename);
//! Read tensor from stream (binary)
    int Read(std::ifstream&);
/*! \brief change upper dimension of tensor, rest zero

The upper structure is resized, For example A(2,3,5) after resize with 7
becomes A(7,0,0) 
 */
    i_number Resize(i_number dd);
/*! \brief Change upper dimension of tensor, preserving available substructure.

The upper structure is resized, For example A(2,3,5) after resize with 7
becomes A(7,...) where substructures A[0](3,5) and A[1](3,5) are preserved, and
new substructures A[2] A[3], ... are created as A[2](0,0). Shrinking is also 
possible, elements that do not fit are deleted.  
 */
    i_number Expand(i_number dd); 
    /*! \brief Change all dimensions of tensor

      The tensor dimensions are changed, e.g. A(2,8,1).FullResize(1,7,5) 
      becomes A(1,7,5), all elements are created UNASSIGNED. 
    */
    i_number FullResize(i_number dd, ...);
      /*! \brief Change all dimensions of tensor and preserve elements
	In addition to  FullResize(i_number dd, ...), all allowed 
	elements are copied, extra created unassigned. 
       */
    i_number FullExpand(i_number *);
    i_number FullExpand(i_number dd, ...);
      //!return rank dimensions in array A(1,6,8).Dim(int *a) gives a={1,6,8} 
    i_number Dim (i_number *); 
    
    int TimesVector(cnumber *yl, cnumber *y, cnumber *x, cnumber shift);
    int TimesVector(cnumber *y, cnumber *x, cnumber shift);


   };

/// \cond
    /*! \class Tensor<1 cnumber>
    \brief 1-dim vectors, base class for higher rank tensors
    */
  template <typename cnumber>
  class Tensor<1,cnumber>
  {
  private:
  
  public: 
    i_number dim;
    cnumber *elem;
    Tensor(); //default constructor 
    Tensor(i_number*); //main constructor
    Tensor(i_number); //create tensor of a given single dimension
    //Tensor(int a) {Tensor(i_number(a));}; //compatibility make from int, creation problem EE(n)???
    Tensor(const char *filename) {Read(filename);}; //reading from file
    //!constractor from itself 
    Tensor(const Tensor<1,cnumber> &); 
    ~Tensor(void); //destructor
    inline operator cnumber* (void) const {return elem;};
        inline i_number size(void) {return dim;}  ///< size function for compatibility 

    Tensor<1,cnumber>& operator=(const Tensor<1,cnumber>& T); 
    template <typename qnumber>
    Tensor<1,cnumber>&  operator=(const Tensor<1,qnumber>& T); 
 //Operators Tensor A+=number 
 template <typename qnumber>
    Tensor<1,cnumber>&  operator *= (const qnumber T); 
 template <typename qnumber>
    Tensor<1,cnumber>&  operator /= (const qnumber T);
 template <typename qnumber> 
    Tensor<1,cnumber>&  operator = (const qnumber T);
 template <typename qnumber>
    Tensor<1,cnumber>&  operator += (const qnumber T);
 template <typename qnumber>
    Tensor<1,cnumber>&  operator -= (const qnumber T);
 
    int Write(const char *filename);
    int Write(std::ofstream &);
    int Read(const char *filename);
   int Read(std::ifstream&);
    //friend std::ostream& operator << (std::ostream&, const Tensor &);

    i_number Resize(i_number); /*! change dimension of tensor*/
    i_number Expand(i_number); /*! change dimension of tensor but do 
				 not destroy*/
    i_number FullExpand(i_number *); /*! change dimension of tensor but do 
				 not destroy*/

    i_number Dim(i_number *); /*return dimension of tensor in array*/
    };


    ///\endcond

//############## Constructors and Destructors #################//

  /*This is a constructor that identifies rank*/
 template <int rank, typename cnumber>  
 Tensor<rank, cnumber>::Tensor() 
  {
     dim=(global_pointer-rank+1)[0];
     //   cout<<"creating tensor rank "<<rank<<" "<<dim<<endl;
    try {elem=new Tensor<rank-1,cnumber>[dim];}
	//catch-throw
	 catch (std::exception& ex)
    {throw ex;}
     }  

  template <int rank, typename cnumber>  
  Tensor<rank, cnumber>::Tensor(i_number *dmn) 
  {
  this->initialize(dmn); //call initialization
   }  

  template <int rank, typename cnumber>  
  Tensor<rank, cnumber>::Tensor(i_number d1, ...):dim(d1) 
  {
    va_list ap;
    va_start(ap,d1); //initialize list
    i_number *dmn=new i_number [rank]; //create an array of remaining dim
    dmn[0]=dim;
    for (int i=1;i<rank;i++) dmn[i]=va_arg(ap,i_number);// read dimensions
    va_end(ap);
	this->initialize(dmn); //call initialization
	delete [] dmn;
     }  
 template <int rank, typename cnumber>  
 Tensor<rank, cnumber>::Tensor(const Tensor<rank, cnumber> &T) 
  {
    dim=T.dim;
	try {elem=new Tensor<rank-1,cnumber>[dim];}
	//catch-throw
	 catch (std::exception& ex)
    {
	FatalError("Exception: constructor Tensor(&Tensor) "<<ex.what());
	}
    for(int i=0;i<dim;i++) elem[i]=T.elem[i];
     }  

  template <int rank, typename cnumber>  
  i_number Tensor<rank, cnumber>::FullResize(i_number dd, ...) 
  {
    va_list ap;
    va_start(ap,dd); //initialize list
    if (dim!=dd) { //only do something if tensor changes
    i_number *dmn=new i_number [rank]; //create an array of remaining dim
    dmn[0]=dd;
    for (int i=1;i<rank;i++) dmn[i]=va_arg(ap,i_number);// read dimensions
    va_end(ap);
    global_pointer=&dmn[rank-1];
    //this is first call and we store global pointer to last dimension
    dim=(global_pointer-rank+1)[0];
    if (elem) delete [] elem; //delete old element
    elem=new Tensor<rank-1,cnumber>[dim];
    delete [] dmn;
   global_pointer=(global_zero+31);
    }

    return dim;
     }  

 template <int rank, typename cnumber>  
   i_number Tensor<rank, cnumber>::FullExpand(i_number *dmn) 
  {
    //very inefficient
 
    if (dmn[0]!=dim) {
    Tensor<rank-1,cnumber> *elemtmp=new Tensor<rank-1,cnumber>[dmn[0]];
   i_number ldim=Min(dmn[0],dim); 
    for (int i=0;i<ldim;i++) Swap(elem[i],elemtmp[i]);
    Swap(elem,elemtmp);    
    if (elemtmp) delete [] elemtmp; //delete old element
    dim=dmn[0];
    }
    for (int i=0;i<dim;i++) elem[i].FullExpand((dmn+1));
    return dmn[0];
  }  



 template <int rank, typename cnumber>  
  i_number Tensor<rank, cnumber>::FullExpand(i_number dd, ...) 
  {
    va_list ap;
    va_start(ap,dd); //initialize list
    i_number *dmn=new i_number [rank]; //create an array of remaining dim
    dmn[0]=dd;
    for (int i=1;i<rank;i++) dmn[i]=va_arg(ap,i_number);// read dimensions
    va_end(ap);
    FullExpand(dmn);
    if (dmn) delete [] dmn;

    return dim;
     }  

 template <int rank, typename cnumber>  
  i_number Tensor<rank, cnumber>::Resize(i_number dd) 
  {
    if (dim!=dd) { //only do something if tensor changes
    i_number *dmn=new i_number [rank]; //create an array of remaining dim
    dmn[0]=dd;
    for (int i=1;i<rank;i++) dmn[i]=0;// set dimensions of sub parts to zero
    global_pointer=&dmn[rank-1];
    //this is first call and we store global pointer to last dimension
    dim=(global_pointer-rank+1)[0];
    if (elem) delete [] elem; //delete old element
    elem=new Tensor<rank-1,cnumber>[dim];
    delete [] dmn;
   global_pointer=(global_zero+31);
    }

    return dim;
     }  

   template <int rank, typename cnumber>  
  i_number Tensor<rank, cnumber>::Expand(i_number dd) 
  {
    if (dim!=dd) { //only do something if tensor changes
    i_number *dmn=new i_number [rank]; //create an array of remaining dim
    dmn[0]=dd;
    for (int i=1;i<rank;i++) dmn[i]=0;// set dimensions of sub parts to zero
    global_pointer=&dmn[rank-1];
    //this is first call and we store global pointer to last dimension
   
    Tensor<rank-1,cnumber> *elemtmp=new Tensor<rank-1,cnumber>[dd];
    i_number ldim=Min(dd,dim); 
    for (int i=0;i<ldim;i++) Swap(elem[i],elemtmp[i]);
    Swap(elem,elemtmp);
    dim=(global_pointer-rank+1)[0]; //overwrite dim
    if (elemtmp) delete [] elemtmp; //delete old element
    delete [] dmn; //delete array matrix
    global_pointer=(global_zero+31); //restore zero pointer
    }

    return dim;
     }  


 template <int rank, typename cnumber>  
   i_number Tensor<rank, cnumber>::Dim(i_number *dmn) 
  {
    if (dim>0) elem[0].Dim(dmn+1); 
    else for (int i=1;i<rank;i++) dmn[i]=0; //get everything to zero
    return dmn[0]=dim;
  }  


 template <int rank, typename cnumber>  
  Tensor<rank, cnumber>::~Tensor(void) 
  {
   delete [] elem;
   elem=NULL;
     }  


 template <int rank, typename cnumber>
 void Swap(Tensor<rank, cnumber>& A, Tensor<rank, cnumber>& B )
 {
   Swap(A.dim, B.dim); 
   Swap(A.elem,B.elem); //call of swap pointers
  } 

/*
  template <int rank, typename cnumber>
  Tensor<rank-1, cnumber>& Tensor<rank, cnumber>::operator [] (i_number n)
  {return elem[n];}
*/
 template <int rank, typename cnumber>
  Tensor<rank, cnumber>& Tensor<rank, cnumber>::operator=(const Tensor<rank, cnumber>& T)
 {
   Resize(T.dim);
      for (int i=0;i<dim;i++) elem[i]=T.elem[i]; //recursive call should be ok with dif types

    return *this;
  } 

  template <int rank, typename cnumber> template<typename qnumber>
  Tensor<rank, cnumber>& Tensor<rank, cnumber>::operator=(const Tensor<rank, qnumber>& T)
  {
    
    Resize(T.dim);
      for (int i=0;i<dim;i++) elem[i]=T.elem[i]; //recursive call should be ok with dif types
    
    return *this;
  }

  template <int rank, typename cnumber>
  int Tensor<rank, cnumber>::Write(const char *filename) {
  std::ofstream outf;
  outf.open(filename, std::ios::binary);
  outf.write((char*)(&dim),sizeof(i_number)); //we save rank;
  for (i_number i=0;i<dim;i++) elem[i].Write(outf);
  outf.close();
 return 1;
  }

 template <int rank, typename cnumber>
 int Write(Tensor<rank, cnumber> &T, const char *filename) {
  std::ofstream outf;
  outf.open(filename, std::ios::binary);
  outf.write((char*)(&T.dim),sizeof(i_number)); //we save rank;
  for (i_number i=0;i<T.dim;i++) Write(T.elem[i],outf);
  outf.close();
 return 1;
  }
  
  template <int rank, typename cnumber>
  int Tensor<rank, cnumber>::Write(std::ofstream &outf) {
   outf.write((char*)(&dim),sizeof(i_number)); //we save rank;
    for (i_number i=0;i<dim;i++) elem[i].Write(outf);
    return 1;
  }

   template <int rank, typename cnumber>
   int Write(Tensor<rank, cnumber> &T, std::ofstream &outf) {
   outf.write((char*)(&T.dim),sizeof(i_number)); //we save rank;
   for (i_number i=0;i<T.dim;i++) Write(T.elem[i],outf);
    return 1;
  }

 template <int rank, typename cnumber>
  int Tensor<rank, cnumber>::Read(const char *filename) {
  std::ifstream inf;
  inf.open(filename, std::ios::binary);
  int dd;
  inf.read((char*)(&dd),sizeof(i_number)); //we read rank;
  Resize(dd);
  //we manipulate hare if rank does not mach
  for (i_number i=0;i<dim;i++) elem[i].Read(inf);
  inf.close();
 return 1;
  }
  
template <int rank, typename cnumber>
int Read(Tensor<rank, cnumber> &T, std::ifstream &inf) {
        int dd;
        //Read(dd, inf);  //we read rank;
        inf.read((char*)(&dd),sizeof(int));
        T.Resize(dd);
   //manipulate here to assure dimension
        for (i_number i=0;i<T.dim;i++) Read(T.elem[i], inf);
        return 1;
      }
      
 template <int rank, typename cnumber>
  int Read( Tensor<rank, cnumber> &T, const char *filename) {
  std::ifstream inf;
  inf.open(filename, std::ios::binary);
  if (!inf.good()) FatalError("Unable to open: "<<filename);
  int dd;
 // Read(dd,inf); //we read rank;
   inf.read((char*)(&dd),sizeof(int));
  T.Resize(dd); //we resize vector
  //we manipulate hare if rank does not mach
  for (i_number i=0;i<T.dim;i++) Read(T.elem[i], inf);
  inf.close();
 return 1;
  }

template <int rank, typename cnumber>
  int Tensor<rank, cnumber>::Read(std::ifstream &inf) {
   int dd;
  inf.read((char*)(&dd),sizeof(i_number)); //we read rank;
  Resize(dd);
   //manipulate here to assure dimension
  for (i_number i=0;i<dim;i++) elem[i].Read(inf);
  return 1;
  }



  template<int rank, typename cnumber>
std::ostream& operator << (std::ostream &stream,
			   Tensor<rank, cnumber> &T)
  {
    for (i_number i=0;i<T.dim;i++) stream<<T.elem[i]<<endl;
    return stream;
  }

  template<int rank, typename cnumber>
  std::istream& operator >> (std::istream &stream,
			   Tensor<rank, cnumber> &T)
  {
    for (i_number i=0;i<T.dim;i++) stream>>T.elem[i];
    return stream;
  }


 template <int rank,typename cnumber> template<typename qnumber>
  Tensor<rank, cnumber>& Tensor<rank, cnumber>::operator *=(const qnumber x)
  {
    for (int i=0;i<dim;i++) elem[i]*=cnumber(x); 
    return *this; 
  }
 template <int rank,typename cnumber> template<typename qnumber>
  Tensor<rank, cnumber>& Tensor<rank, cnumber>::operator /=(const qnumber x)
  {
    for (int i=0;i<dim;i++) elem[i]/=cnumber(x); 
    return *this; 
  }
   template <int rank,typename cnumber> template<typename qnumber>
  Tensor<rank, cnumber>& Tensor<rank, cnumber>::operator =(const qnumber x)
  {
    for (int i=0;i<dim;i++) elem[i]=cnumber(x); 
    return *this; 
  }
   template <int rank,typename cnumber> template<typename qnumber>
  Tensor<rank, cnumber>& Tensor<rank, cnumber>::operator +=(const qnumber x)
  {
    for (int i=0;i<dim;i++) elem[i]+=cnumber(x); 
    return *this; 
  }
   template <int rank,typename cnumber> template<typename qnumber>
  Tensor<rank, cnumber>& Tensor<rank, cnumber>::operator -=(const qnumber x)
  {
    for (int i=0;i<dim;i++) elem[i]-=cnumber(x); 
    return *this; 
  }
 
  //end operator Tensor(rank...) += cnumber

  //left-right multiplication for tensor rank 2 is assumed.
  template <int rank, class cnumber>
  int Tensor<rank,cnumber>::TimesVector(cnumber *y, cnumber *x, cnumber shift)
{

  for (int i=0;i<dim;i++) {
    y[i]=-shift*x[i]; //zero vector y
    for (int j=0;j<dim;j++)  y[i]+=elem[i][j]*x[j];
  }

 return dim;
}

  

  //left-right multiplication for tensor rank 2 is assumed No Complex NUMBERS.
  template <int rank, class cnumber>
  int Tensor<rank,cnumber>::TimesVector(cnumber *yl, cnumber *y, cnumber *x, cnumber shift)
{

 for (int i=0;i<dim;i++) y[i]=-shift*x[i]; //zero vector y
 for (int i=0;i<dim;i++) yl[i]=-(shift*x[i]); //zero vector y -(shift|x[i]); 
 for (int i=0;i<dim;i++) for (int j=0;j<dim;j++)  
   {                 
		  y[i]+=elem[i][j]*x[j];
                  yl[i]+=(elem[j][i]*x[j]); //times conjugate (elem[j][i]|x[j]);
 }

 return dim;
}
  /* this is for hermitian
  template <int rank>
  int Tensor<rank,Complex>::TimesVector(Complex *yl, Complex *y, Complex *x, Complex shift)
{

 for (int i=0;i<dim;i++) y[i]=-shift*x[i]; //zero vector y
 for (int i=0;i<dim;i++) yl[i]=-(shift|x[i]); //zero vector y
 for (int i=0;i<dim;i++) for (int j=0;j<dim;j++)  
   {                 
		  y[i]+=elem[i][j]*x[j];
                  yl[i]+=(elem[j][i]|x[j]); //times conjugate
 }

 return dim;
}
  */ //end comment



  /*THIS PART IS FOR RANK 1 */
///\cond 
template <typename cnumber>  
  Tensor<1, cnumber>::Tensor() 
  {
     dim=(global_pointer)[0];
    try {elem=new cnumber [dim];}
	//throw back
	 catch (std::exception& ex) {throw ex;}
     }  

template <typename cnumber>
Tensor<1, cnumber>::Tensor(i_number *dmn) 
  {
    dim=dmn[0];
    try {elem=new cnumber [dim];}
	catch (std::exception& e)
    {
	cerr<<"Exception: constructor Tensor<1>(&Tensor(";
	cerr<<dim<<"))"<<endl;
    FatalError(e.what());
    }
	
  }

template <typename cnumber>
Tensor<1, cnumber>::Tensor(const Tensor<1,cnumber> &T) 
  {
    dim=T.dim;
    elem=new cnumber [dim];
    for (int i=0;i<dim;i++) elem[i]=T.elem[i];
  }

template <typename cnumber>
Tensor<1, cnumber>::Tensor(i_number dmn):dim(dmn) 
  {
     try { elem=new cnumber [dim];}
	  catch (std::exception& e)
    {
	cerr<<"Exception: constructor Tensor<1>(";
	cerr<<dim<<")"<<endl;
    FatalError(e.what());
    }
  }

template <typename cnumber>
Tensor<1, cnumber>::~Tensor(void) 
  {
    delete [] elem;
  }

template <typename cnumber>
i_number Tensor<1, cnumber>::Resize(i_number dmn) 
  {
    if (dmn!=dim) { //only work if dimension are not equal
      dim=dmn;
      if(elem) delete [] elem;
      try {elem=new cnumber [dim];} catch (std::exception& e)
      {
        cerr<<"Exception: function Tensor<1>.Resize(";
        cerr<<dim<<")"<<endl;
        FatalError(e.what());
      }
    }
    return dim;
  }

template <typename cnumber>
i_number Tensor<1, cnumber>::Expand(i_number dmn) 
  {
    if (dmn!=dim) { //only work if dimension are not equal
      //dim=dmn;
      cnumber *elemtmp=new cnumber [dmn];
      i_number ldim=Min(dim,dmn); 
      for (int i=0;i<ldim;i++) elemtmp[i]=elem[i];
      if(elem) delete [] elem;
      elem=elemtmp;
      dim=dmn;
}
    return dim;
  }

template <typename cnumber>
i_number Tensor<1, cnumber>::Dim(i_number *dmna) 
  {
    
    return dmna[0]=dim; //return trivial dimension
  }


template <typename cnumber>
i_number Tensor<1, cnumber>::FullExpand(i_number *dmna) 
  {
    
    return Expand(dmna[0]); //expand from array
  }


  template <typename cnumber> 
  Tensor<1, cnumber>& Tensor<1, cnumber>::operator=(const Tensor<1, cnumber>& T)
  {
       //if (this != &T) 
    Resize(T.dim);
          for (int i=0;i<dim;i++) elem[i]=T.elem[i]; //recursive call 
    return *this; //multiple assignment possible
  }

  template <typename cnumber> template<typename qnumber>
  Tensor<1, cnumber>& Tensor<1, cnumber>::operator=(const Tensor<1, qnumber>& T)
  {
     //if (this != &T) 
      Resize(T.dim);
      for (int i=0;i<dim;i++) elem[i]=(cnumber)T.elem[i]; //recursive call 
    
    return *this; //multiple assignment possible
  }

template <typename cnumber>
  int Tensor<1, cnumber>::Write(const char *filename) {
  std::ofstream outf;
  outf.open(filename, std::ios::binary);
  outf.write((char*)(&dim),sizeof(i_number)); //we save rank;
  outf.write((char*)(elem), dim*sizeof(cnumber));
 outf.close();
 return 1;
  }


template <typename cnumber>
int Write(Tensor<1, cnumber> &T, const char *filename) {
  std::ofstream outf;
  outf.open(filename, std::ios::binary);
  outf.write((char*)(&T.dim),sizeof(i_number)); //we save rank;
  for (i_number i=0;i<T.dim;i++) Write(T.elem[i],outf);
 outf.close();
 return 1;
  }

  
  template <typename cnumber>
  int Tensor<1, cnumber>::Write(std::ofstream &outf) {
   outf.write((char*)(&dim),sizeof(i_number)); //we save rank;
    outf.write((char*)(elem), dim*sizeof(cnumber));
  return 1;
  }

    template <typename cnumber>
    int Write(Tensor<1, cnumber> &T, std::ofstream &outf) {
      outf.write((char*)(&T.dim),sizeof(i_number)); //we save rank;
      for (i_number i=0;i<T.dim;i++) Write(T.elem[i],outf);
  return 1;
  }

template <typename cnumber>
  int Tensor<1, cnumber>::Read(const char *filename) {
  std::ifstream inf;
  inf.open(filename, std::ios::binary);
  i_number rdim;
  inf.read((char*)(&rdim),sizeof(i_number)); //we read rank;
  Resize(rdim); //resize tensor before read
  inf.read((char*)(elem), dim*sizeof(cnumber));
  inf.close();
 return 1;
  }
  
template <typename cnumber>
int Read(Tensor<1, cnumber> &T, const char *filename) {
  std::ifstream inf;
  inf.open(filename, std::ios::binary);
  i_number rdim;
  inf.read((char*)(&rdim),sizeof(i_number)); //we read rank;
  T.Resize(rdim); //resize tensor before read
  for (i_number i=0;i<T.dim;i++) Read(T.elem[i], inf);
  inf.close();
 return 1;
  }


  template <typename cnumber>
  int Tensor<1, cnumber>::Read(std::ifstream &inf) {
   i_number rdim;
   inf.read((char*)(&rdim),sizeof(i_number)); //we read rank;
   //need to change dimensions appropriately
   Resize(rdim); //resize tensor before read
   inf.read((char*)(elem), dim*sizeof(cnumber));
  return 1;
  }

 template <typename cnumber>
 int Read(Tensor<1, cnumber> &T, std::ifstream &inf) {
   i_number rdim;
   inf.read((char*)(&rdim),sizeof(i_number)); //we read rank;
   //need to change dimensions appropriately
   T.Resize(rdim); //resize tensor before read
   //   inf.read((char*)(elem), dim*sizeof(cnumber));
   for (i_number i=0;i<T.dim;i++) Read(T.elem[i], inf);
  return 1;
  }

template<typename cnumber>
std::ostream& operator << (std::ostream &stream,
			   Tensor<1, cnumber> &T)
  {
    for (i_number i=0;i<T.dim;i++) stream<<T.elem[i]<<" ";
    return stream;
  }

template<typename cnumber>
std::istream& operator >> (std::istream &stream,
			   Tensor<1, cnumber> &T)
  {
    for (i_number i=0;i<T.dim;i++) stream>>T.elem[i];
    return stream;
  }

  //*************************Operators************************
  template <typename cnumber,typename qnumber>  
  cnumber operator * (const Tensor<1,cnumber> &a, const Tensor<1,qnumber> &b)
{
        
  //if (a.m!=b.n)
  //              FatalError (" 'Matrix*Matrix':\nwrong dimentions");
        cnumber ret=0;
	//cnumber* ap=a.elem;
	//qnumber* bp=b.elem;
	//cnumber* aendp=a.elem+a.dim;
	//do ret+=(*(ap++))*(*(bp++)); while (ap!=aendp); //10.92
	//	for (int i=0;i<a.dim;i++)  ret+=(*(ap++))*(*(bp++)); //3.6 
	for (int i=0;i<a.dim;i++)  ret+=a.elem[i]*b.elem[i];  //10.85(5) per 30M  
        return ret;
}







 template <typename cnumber> template<typename qnumber>
  Tensor<1, cnumber>& Tensor<1, cnumber>::operator *=(const qnumber x)
  {
    for (int i=0;i<dim;i++) elem[i]*=cnumber(x); 
    return *this; 
  }
 template <typename cnumber> template<typename qnumber>
  Tensor<1, cnumber>& Tensor<1, cnumber>::operator /=(const qnumber x)
  {
    for (int i=0;i<dim;i++) elem[i]/=cnumber(x); 
    return *this; 
  }
   template <typename cnumber> template<typename qnumber>
  Tensor<1, cnumber>& Tensor<1, cnumber>::operator =(const qnumber x)
  {
    for (int i=0;i<dim;i++) elem[i]=cnumber(x); 
    return *this; 
  }
   template <typename cnumber> template<typename qnumber>
  Tensor<1, cnumber>& Tensor<1, cnumber>::operator +=(const qnumber x)
  {
    for (int i=0;i<dim;i++) elem[i]+=cnumber(x); 
    return *this; 
  }
   template <typename cnumber> template<typename qnumber>
  Tensor<1, cnumber>& Tensor<1, cnumber>::operator -=(const qnumber x)
  {
    for (int i=0;i<dim;i++) elem[i]-=cnumber(x); 
    return *this; 
  }
 
///\endcond

} //end namespace
#endif
	  

