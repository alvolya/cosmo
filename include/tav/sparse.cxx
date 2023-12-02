/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/** 
\author Alexander Volya  <http://www.volya.net>
\file sparse.cxx
\brief Sparse matrix class, indexed rows of non-zero elements 
*/
/*
09/05/01 new class of symmetric matreces from sparse, it is very similar
11/04/2005 forcefully write element, WriteElement(i,j)
04/25/06 include variable type=0 normal, 1-symmetr
access to element uses swap
04/19/2006 invalid constructor function is fixed g++ v4
07/17/07 include left-right multiplication, no shift
11/13/07 bugfix no delete matrix
01/28/08 An important bug is fixed in Sync(); 
06/11/08 Compatibility with Matrix is removed.  
07/14/08 Times vector templated, no longer pointer any array with [] access
07/3/11 TimesMatrix is included for davidson method
*/
#ifndef __TAV_SPARSE_CXX__
#define __TAV_SPARSE_CXX__
#include <cstdlib>
#include <fstream>
#include <iostream>
//#include "basic.cxx"
//using tav::Abs;
//#include "matrix.cxx"

#define dimension_of_space 0;
#ifndef DPRECISION
#define DPRECISION 1E-8;
#endif

#define FatalError_sp(message)              { std::cerr <<"Error: "<< message <<std::endl<< __FILE__ << ":" << __LINE__ << std::endl; exit(1); }

namespace tav{
    typedef unsigned long long  uint_nbasis;
/** \brief Sparse matrix class
  */
template <class cnumber=double>
class SparseMatrix
{ 
  /* structure consists of following
   - matrix size is n*m,
   - indexing works by rows n
   - nz[i] is the number of nonzero elemets in row i
   - ind[i][j] gives for the row i, the column number of j-th nonzero elemet
   - elem[i][j] gives the numerical value of corresponding to ind[i][j] elemet
  */
  /*  protected:
  cnumber *elem;
  int    *nz;
  int    **ind; */
  public:
  int type; ///<0-normal 1-symmetric
  cnumber **elem; ///<non-zero matrix elements elem[i] points on i-th row
  unsigned    *nz; ///<nz, number of non zero elements in nz[i]
  uint_nbasis    **ind; ///< list of j-coordinates of non-zero elements. 
  cnumber zeroelement;  //this controls zero elemets, allows to use reference
  uint_nbasis zeroi,zeroj;     
  double precision; ///<precision with which matrix elements are rejected
  //upper part for testing only
  uint_nbasis n;                  ///<  n-size of the matrix  
  SparseMatrix(void);     //  create a defult zero matrix //
  SparseMatrix(uint_nbasis);
  template <class rnumber> 
  SparseMatrix(SparseMatrix<rnumber> &); //constract 
  ~SparseMatrix(void);   //delete matrix
  cnumber& operator () (uint_nbasis, uint_nbasis);
  cnumber& WriteElement(uint_nbasis, uint_nbasis); //full access of element non-symmetric
  //  double operator () (int, int) const;
  uint_nbasis Sync (void);
  void WriteLine(uint_nbasis i, unsigned nz, uint_nbasis *ind ,cnumber* elem);
 template <class qnumber> 
  void TimesVector(qnumber& ,qnumber&);
template <class qnumber> 
void TimesVector(qnumber &, qnumber &,qnumber &); //left-right multiplication
  template <class qnumber> 
 void TimesVector(qnumber&,qnumber&, cnumber);
  void TimesMatrix(cnumber**,cnumber**, unsigned) const; //<multiply many vectors at the same time
 void TimesMatrix(cnumber*,cnumber*, unsigned) const; //<multiply many vectors at the same time
 uint_nbasis ReadDiagonal(cnumber *x); //read diagonal
  //  void TimesVector(cnumber*,cnumber*);
  template <class rnumber>
  SparseMatrix operator = (SparseMatrix<rnumber> &);
   uint_nbasis Write(const char *filename);
  uint_nbasis Read(const char *filename);
  uint_nbasis size(){return n;};
  //template <class rnumber> 
  //rnumber& operator [] (int);

  //matrix compatibility
   // SparseMatrix operator = (const Matrix<cnumber> &);
 //SparseMatrix(const Matrix<cnumber> &); //from double_matrix  

};

  //CONSTRUCTORS++++++++++++++++++++++++++++++++++++++++

template <class cnumber>
SparseMatrix<cnumber>::SparseMatrix(void)
// by default create zero matrix of size 4
{
  zeroelement=0.0;
  n=dimension_of_space;
    zeroi=n;
    zeroj=n;
  nz=new unsigned[n];
    ind=new uint_nbasis*[n];
    elem=new cnumber*[n];
  for (uint_nbasis i=0;i<n;i++) {nz[i]=0; ind[i]=NULL; elem[i]=NULL;}
  precision=1E-10;
}

template <class cnumber>
SparseMatrix<cnumber>::SparseMatrix(uint_nbasis v1):n(v1)
  //matrix with defined dimenstion
{   zeroelement=0.0;
    zeroi=n;
    zeroj=n;
precision=1E-10;
  nz=new unsigned[n];
  ind=new uint_nbasis*[n];
   elem=new cnumber*[n];
  for (uint_nbasis i=0;i<n;i++) {nz[i]=0; ind[i]=NULL; elem[i]=NULL;}
}

  template <class cnumber> template <class rnumber>
SparseMatrix<cnumber>::SparseMatrix(SparseMatrix<rnumber> &mat):n(mat.n)
  //create matrix from other sparse matrix
{ 
  if (std::abs(mat.zeroelement)>precision) mat.Sync();
  zeroelement=0.0;
    zeroi=n;
    zeroj=n;
precision=1E-10;
  nz=new unsigned[n];
  ind=new uint_nbasis*[n];
  elem=new cnumber*[n];
 // int num=0;
  for (uint_nbasis i=0;i<n;i++) {
    nz[i]=mat.nz[i];        //copy number of non-zeros for each row
  ind[i]=new uint_nbasis [nz[i]];
  elem[i]=new cnumber [nz[i]];
  for (unsigned j=0;j<nz[i];j++) {
    ind[i][j]=mat.ind[i][j]; //copy indexing
    elem[i][j]=cnumber(mat.elem[i][j]); //copy elemets
  }  
  }
}
template <class cnumber>
SparseMatrix<cnumber>::~SparseMatrix(void)
{ 
        for (uint_nbasis i=0;i<n;i++) if (ind[i]) {
         delete []ind[i];
         delete []elem[i];
         ind[i]=NULL;
         elem[i]=NULL;
	}
        /*
        if (nz)
        {
                delete []nz;
                nz=NULL;
        }
        */
        delete [] nz;
        delete [] ind;
        delete [] elem;
}

/*Accessing elemets   */
template <class cnumber>
cnumber& SparseMatrix<cnumber>::WriteElement(uint_nbasis ii,uint_nbasis jj)
  /* this will work for x=A(i,j), getting matrix element */
{       
  if (std::abs(zeroelement)>precision) Sync();
        if ((ii>=n)||(jj>=n)) FatalError_sp ("SparseMatrix::operator(int,int)\nBad arguments.\n");
	uint_nbasis i,j;
        //if (ii<=jj) {i=ii;j=jj;}  else {i=jj; j=ii;};  
   for (unsigned q=0;q<nz[ii];q++)
     { 
   if (ind[ii][q]==jj) return elem[ii][q]; }
   /*this matrix elemet is not in the list*/
   zeroi=ii;
   zeroj=jj;
   return zeroelement;
}

template <class cnumber>
cnumber& SparseMatrix<cnumber>::operator () (uint_nbasis ii,uint_nbasis jj)
  /* this will work for x=A(i,j), getting matrix element */
{       
  if (std::abs(zeroelement)>precision) Sync();
        if ((ii>=n)||(jj>=n)||(jj<0)||(ii<0)) FatalError_sp ("SparseMatrix::operator(int,int)\nBad arguments.\n");
	//	uint_nbasis i,j;
	//        if (ii<=jj) {i=ii;j=jj;}  else {i=jj; j=ii;};  
	if ((type)&&(ii>jj)) std::swap(ii,jj);
   for (unsigned q=0;q<nz[ii];q++)
     { 
   if (ind[ii][q]==jj) return elem[ii][q]; }
   /*this matrix elemet is not in the list*/
   zeroi=ii;
   zeroj=jj;
   return zeroelement;
}


template <class cnumber>
void  SparseMatrix<cnumber>::WriteLine(uint_nbasis line,unsigned n, uint_nbasis *in, cnumber *el) {
    if (std::abs(zeroelement)>precision) Sync();
    nz[line]=n;
    if (ind[line]) delete []ind[line];
    if (elem[line]) delete []elem[line]; 
  ind[line]=new uint_nbasis [n];
  elem[line]=new cnumber [n];
  for (uint_nbasis i=0;i<n;i++) {ind[line][i]=in[i]; elem[line][i]=el[i];} 
}

template <class cnumber> template <class qnumber> 
void  SparseMatrix<cnumber>::TimesVector(qnumber & c,qnumber &a) {
  // c= matrix*a
    if (std::abs(zeroelement)>precision) Sync();
	for (uint_nbasis i=0;i<n;i++) c[i]=0.0;
        for (uint_nbasis i=0;i<n;i++) {
          
                for (unsigned k=0;k<nz[i];k++) {
		  c[i]+=elem[i][k]*a[ind[i][k]];
                  if (type) if(i!=ind[i][k]) c[ind[i][k]]+=elem[i][k]*a[i];}
	}
}


template <class cnumber> template <class qnumber> 
void  SparseMatrix<cnumber>::TimesVector(qnumber & c1, qnumber & c, qnumber &a) {
  // c= matrix*a
    if (std::abs(zeroelement)>precision) Sync();
	for (uint_nbasis i=0;i<n;i++) c1[i]=c[i]=0.0;
        for (uint_nbasis i=0;i<n;i++) {
          
                for (unsigned k=0;k<nz[i];k++) {
		  c[i]+=elem[i][k]*a[ind[i][k]];
                  if (type) if(i!=ind[i][k]) c[ind[i][k]]+=elem[i][k]*a[i];}

		//multiplication by transpose
		for (unsigned k=0;k<nz[i];k++) {
		  c1[ind[i][k]]+=elem[i][k]*a[i];
                  if (type) if(i!=ind[i][k])  c1[i]+=elem[i][k]*a[ind[i][k]];}

	}
}

template <class cnumber> template <class qnumber> 
void  SparseMatrix<cnumber>::TimesVector(qnumber &c,qnumber &a, cnumber shift) {
  // c= matrix*a
  //cerr<<"Multiplication with shift : "<<shift<<endl;
    if (std::abs(zeroelement)>precision) Sync();
	for (uint_nbasis i=0;i<n;i++) c[i]=-shift*a[i];
        for (uint_nbasis i=0;i<n;i++) {
          
                for (unsigned k=0;k<nz[i];k++) {
		  c[i]+=elem[i][k]*a[ind[i][k]];
                  if (type) if(i!=ind[i][k]) c[ind[i][k]]+=elem[i][k]*a[i];}
	}
}
//note for compatibility with davidson this must be constant
template <class cnumber>  
void  SparseMatrix<cnumber>::TimesMatrix(cnumber **yy, cnumber **xx, unsigned m) const {
  // c= matrix*a
  //cerr<<"Multiplication with shift : "<<shift<<endl;
    if (std::abs(zeroelement)>precision) FatalError_sp("no sync done before SparseMatrix::TimesMatrix");

 
#pragma omp parallel for schedule(dynamic, 1)
    	for(unsigned im=0;im<m;im++) {//loop over m vectors
        for (uint_nbasis i=0;i<n;i++) yy[im][i]=0.0; //zero vectors y
        for (uint_nbasis i=0;i<n;i++) { //go over all elements
                for (unsigned k=0;k<nz[i];k++) {
			
				yy[im][i]+=elem[i][k]*xx[im][ind[i][k]];
                  if (type) if(i!=ind[i][k]) yy[im][ind[i][k]]+=elem[i][k]*xx[im][i];
				  } //end looping over vectors
				  } //end looping over non-zero elements 
	} //end looping over all elements
} //end routine
template <class cnumber>  
void  SparseMatrix<cnumber>::TimesMatrix(cnumber *y, cnumber *x, unsigned m) const {
  // c= matrix*a
  //cerr<<"Multiplication with shift : "<<shift<<endl;
    if (std::abs(zeroelement)>precision) FatalError_sp("no sync done before SparseMatrix::TimesMatrix");
	cnumber **yy=new cnumber*[m]; //pointers to vectors
    cnumber **xx=new cnumber*[m]; //pointers to out vectors
 for(unsigned im=0;im<m;im++) {//loop over m vectors
            yy[im]=y+n*im; //shift pointer so yy[i] is a column vector
            xx[im]=x+n*im;}
    TimesMatrix(yy, xx, m);
	 delete [] yy;
 delete [] xx;
} //end routine


template <class cnumber> 
uint_nbasis  SparseMatrix<cnumber>::ReadDiagonal(cnumber *x) {
  // c= matrix*a
  //cerr<<"Multiplication with shift : "<<shift<<endl;
    if (std::abs(zeroelement)>precision) Sync();
        for (uint_nbasis i=0;i<n;i++) {          
		if (ind[i][0]==i) x[i]=elem[i][0]; //this is diagonal (we use upper-right)
        else x[i]=0.0;  			
	}
	return n;
}

template <class cnumber>
uint_nbasis SparseMatrix<cnumber>::Sync(void) {
  if (std::abs(zeroelement)<precision) return 0; //should not call sync

  uint_nbasis *tqr= new uint_nbasis [nz[zeroi]+1];
  cnumber *th=new cnumber [nz[zeroi]+1];
  uint_nbasis i=0;
  for (;i<nz[zeroi];i++) if (ind[zeroi][i]<zeroj) {
    tqr[i]=ind[zeroi][i];
    th[i]=elem[zeroi][i];
  }
  else break;
    tqr[i]=zeroj;
    th[i]=zeroelement;
    for (unsigned j=i;j<nz[zeroi];j++) {
    tqr[j+1]=ind[zeroi][j];
    th[j+1]=elem[zeroi][j];}
    /*write back elements*/
    if (elem[zeroi]) delete [](elem[zeroi]);
    if (ind[zeroi]) delete [](ind[zeroi]);
     nz[zeroi]++;
     elem[zeroi]=th;
     ind[zeroi]=tqr;
     //cerr<<"Call Sync () "<<zeroi<<" "<<zeroj<<" "<<zeroelement<<endl;
     /*prepare to exit */
    zeroi=n;
    zeroj=n;
     zeroelement=0.0;   
     
     //     cout<<"sync called  "<<endl;
  return 0;
}
template <class cnumber> template <class rnumber>
SparseMatrix<cnumber> SparseMatrix<cnumber>::operator = (SparseMatrix<rnumber> &a)
{if (std::abs(a.zeroelement)>a.precision) a.Sync();
 if ((a.n!=n)) FatalError_sp("'SparseMatrix=SparseMatrix':\nmatrices has diffirent dimentions");
  zeroelement=0.0;
    zeroi=n;
    zeroj=n;
  precision=a.precision;
  for (uint_nbasis i=0;i<n;i++) {
    nz[i]=a.nz[i];
  if (ind[i]) delete []ind[i];
  if (elem[i]) delete []elem[i];
  ind[i] = new uint_nbasis [nz[i]];
  elem[i]= new cnumber [nz[i]];
  for (unsigned j=0;j<nz[i];j++) {
    ind[i][j]=a.ind[i][j];
    elem[i][j]=cnumber(a.elem[i][j]); }}
        return *this;
}

// ################## Operators '<<' and '>>' ###################
template <class cnumber>
std::ostream& operator << (std::ostream &stream,  
SparseMatrix<cnumber> &a)
{       if (std::abs(a.zeroelement)>a.precision) a.Sync();
 for (uint_nbasis i=0;i<a.n;i++) {
   // stream<<a.nz[i]<<'\n';
   for (unsigned j=0;j<a.nz[i];j++) {
stream<<i<<" "<<a.ind[i][j]<<" "<<a.elem[i][j];
 stream<<'\n'; }}
       return stream;
 }

template <class cnumber>
uint_nbasis SparseMatrix<cnumber>::Write(const char *filename) {
  std::ofstream outf;
  outf.open(filename, std::ios::binary);
 if (std::abs(zeroelement)>precision) Sync();
 outf.write((char*)(&n),sizeof(uint_nbasis)); //we have saved n
 //outf.write((char*)(nz),n*sizeof(int)); //we saved nz
 for(uint_nbasis i=0;i<n;i++) {
   outf.write((char*)(&(nz[i])),sizeof(unsigned)); //we saved nz
   outf.write((char*)(ind[i]),nz[i]*sizeof(uint_nbasis)); //saved array of ind
   outf.write((char*)(elem[i]),nz[i]*sizeof(cnumber)); //saved array of elem
 }
 outf.close();
 return 1;
}

template <class cnumber>
uint_nbasis SparseMatrix<cnumber>::Read(const char *filename) {
  std::ifstream outf;
  outf.open(filename, std::ios::binary);
if (!outf)   {               // if the file does not exist
    std::cerr<<"File "<<filename<<" is not found"<<std::endl; // return error message
  return 0;
}
//Delete old stuff if any
       for (uint_nbasis i=0;i<n;i++) if (ind[i]) {
         delete []ind[i];
         delete []elem[i];
         ind[i]=NULL;
         elem[i]=NULL;
	}
        for (uint_nbasis i=0;i<n;i++) if (elem[i]) {
         delete []elem[i];
         elem[i]=NULL;
	}
        if (nz)
        {
                delete []nz;
                nz=NULL;
        }
 //Read now new things
 outf.read((char*)(&n),sizeof(uint_nbasis)); //got n
  zeroelement=0.0;
    zeroi=n;
    zeroj=n;
  nz=new unsigned[n];
  // outf.read((char*)(nz),n*sizeof(int)); //we got nz
  ind=new uint_nbasis*[n];
  elem=new cnumber*[n];
 

 for(uint_nbasis i=0;i<n;i++) {
   outf.read((char*)(&(nz[i])),sizeof(unsigned));
   ind[i]=new uint_nbasis [nz[i]];
   elem[i]=new cnumber [nz[i]];
   outf.read((char*)(ind[i]),nz[i]*sizeof(uint_nbasis)); //got ind
   outf.read((char*)(elem[i]),nz[i]*sizeof(cnumber)); //got elem
 }
 outf.close();
 return 1;
}

  //matrix and sparse matrix interface
/*
template <class cnumber>
Matrix<cnumber>  operator * 
(SparseMatrix<cnumber> a, Matrix<cnumber> &vec) {
  // c= matrix*a
  if (a.n!=vec.n) FatalError_sp("SparseMatrix * Matrix : wrong dimentions");
    if (std::abs(a.zeroelement)>a.precision) a.Sync();
    Matrix<cnumber> c(a.n, vec.m);
    Zero(c);
    for (int j=0;j<vec.m;j++)
        for (int i=0;i<a.n;i++) {
	  //this is line multiplication
          
                for (int k=0;k<a.nz[i];k++) {
		  c(i,j)+=a.elem[i][k]*vec(a.ind[i][k],j);
		  if (a.type) if(i!=a.ind[i][k]) c(a.ind[i][k],j)+=a.elem[i][k]*vec(i,j);}
	}
    return c;
}

template <class cnumber>
Matrix<cnumber>  TransMultiply 
(SparseMatrix<cnumber> a, Matrix<cnumber> &vec) {
 // c= matrix*a
  if (a.n!=vec.n) FatalError("SparseMatrix * Matrix : wrong dimentions");
    if (std::abs(a.zeroelement)>a.precision) a.Sync();
    Matrix<cnumber> c(a.n, vec.m);
    Zero(c);
    for (int j=0;j<vec.m;j++)
        for (int i=0;i<a.n;i++) {
	  //this is line multiplication
          
                for (int k=0;k<a.nz[i];k++) {
		  c(i,j)+=a.elem[i][k]*vec(a.ind[i][k],j);
                  if (a.type) if(i!=a.ind[i][k]) c(a.ind[i][k],j)+=a.elem[i][k]*vec(i,j);}
	}
    return c;
}

template<class cnumber>
SparseMatrix<cnumber>::SparseMatrix(const Matrix<cnumber> &mat):n(mat.n)
  //create matrix from other matrix
{ if (mat.n!=mat.m)  ErrorMsg ("Error in SparseMatrix(Double_Matrix)\nDouble matrix is not square\n");
  zeroelement=0.0;
  type=1;
  zeroi=-1;
  zeroj=-1;
  precision=DPRECISION;
  nz=new int[n];
  ind=new int*[n];
  elem=new double*[n];
  int num=0;
  int count;
  for (int i=0;i<n;i++) {
  nz[i]=0; 
  for (int j=0;j<=i;j++) {
    if (fabs(mat(i,j))>precision) nz[i]++;
    if (fabs(mat(i,j)-mat(j,i))>precision)  ErrorMsg ("Error in SparseMatrix(Double_Matrix): Double_Matrix is not symmetric.\n");
}
  ind[i]=new int [nz[i]];
  elem[i]=new double [nz[i]];
  count=0;
  for (int j=0;j<=i;j++) 
    {if (fabs(mat(i,j))>precision) {
   ind[i][count]=j; //copy indexing
   elem[i][count]=mat(i,j);
   count++;
    }
    }  };
}


template<class cnumber>
SparseMatrix<cnumber> SparseMatrix<cnumber>::operator = (const Matrix<cnumber> &a)
{ if ((a.m!=n) || (a.n!=n)) ErrorMsg("'SparseMatrix=Double_Matrix':\nmatrix is not square"); 
  zeroelement=0.0;
  zeroi=-1;
  zeroj=-1;
  precision=DPRECISION;
  int count;
  for (int i=0;i<n;i++) {
  nz[i]=0; 
  type =1; 
  for (int j=0;j<=i;j++){ 
    if (std::abs(a(i,j))>precision) nz[i]++;
    if (std::abs(a(i,j)-a(j,i))>precision) ErrorMsg("Error in 'SparseMatrix=Double_Matrix':\nDouble matrix is not symmetric");
  }
  if (ind[i]) delete []ind[i];
  if (elem[i]) delete []elem[i]; 
  ind[i]=new int [nz[i]];
  elem[i]=new cnumber [nz[i]];
  count=0;
  for (int j=0;j<=i;j++) 
    {if (std::abs(a(i,j))>precision) {
   ind[i][count]=j; //copy indexing
   elem[i][count]=a(i,j);
   count++;
    }
  }  }
  return *this;
}



template <class cnumber>  
Matrix<cnumber> PreconditionerSolve (SparseMatrix<cnumber> &a, const Matrix<cnumber> &b)
{
  // if (b.num==1) return Trans(a)*b(0);
        if (a.n!=b.n)
                FatalError (" 'Preconditioner solve':\nwrong dimentions");
        Matrix<cnumber> c(a.n, 1);
              for (int i=0;i<a.n;i++)
		c.elem[i]=b.elem[i]/a(i,i);
        return c;
}

template <class cnumber>  
Matrix<cnumber> PreconditionerTransSolve (SparseMatrix<cnumber> &a, const Matrix<cnumber> &b)
{
  // if (b.num==1) return Trans(a)*b(0);
        if (a.n!=b.n)
                FatalError (" 'Preconditioner solve':\nwrong dimentions");
        Matrix<cnumber> c(a.n, 1);
              for (int i=0;i<a.n;i++)
		c.elem[i]=b.elem[i]/a(i,i);
        return c;
}


  //THIS FUNCTION IS USED TO CAST MATRICES
  template <class rnumber, class cnumber> 
  int cast(SparseMatrix<rnumber> &b, SparseMatrix<cnumber> & a)
{if (std::abs(a.zeroelement)>a.precision) a.Sync();
  // SparseMatrix <rnumber> b(a.n);
//Delete old stuff if any
       for (int i=0;i<b.n;i++) if (b.ind[i]) {
         delete []b.ind[i];
         delete []b.elem[i];
         b.ind[i]=NULL;
         b.elem[i]=NULL;
	}
        for (int i=0;i<b.n;i++) if (b.elem[i]) {
         delete []b.elem[i];
         b.elem[i]=NULL;
	}
        if (b.nz)
        {
                delete []b.nz;
                b.nz=NULL;
        }
	if (b.elem)
        {
                delete []b.elem;
                b.elem=NULL;
        }
	if (b.ind)
        {
                delete []b.ind;
                b.ind=NULL;
        }
	

	b.n=a.n;
	b.nz=new int[b.n];
	b.ind=new int*[b.n];
	b.elem=new rnumber*[b.n];
  b.zeroelement=0.0;
  b.zeroi=-1;
  b.zeroj=-1;
  b.precision=a.precision;
  for (int i=0;i<b.n;i++) {
    b.nz[i]=a.nz[i];
    // if (b.ind[i]) delete []b.ind[i];
    //if (b.elem[i]) delete []b.elem[i];
  b.ind[i] = new int [b.nz[i]];
  b.elem[i]= new rnumber [b.nz[i]];
  for (int j=0;j<b.nz[i];j++) {
    b.ind[i][j]=a.ind[i][j];
    b.elem[i][j]=rnumber(a.elem[i][j]); }}
        return 1;
}
*/

}


#endif // __TAV_SPARSE_CXX__ //   

