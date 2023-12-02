/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/** \defgroup gp_mbo Many-body operators 
*/
/*! 
\author Alexander Volya  <http://www.volya.net>
\file PPMatrix.cxx 
\brief Definition of Many-Body virtual operator-matrix class 
\ingroup gp_mbo
*/
#ifndef __PPMatrix__CXX__
#define __PPMatrix__CXX__
//#include <omp.h> //only pragmas are used
#include <av/statusclass.cxx>
#include "cioperators.cxx"
#include <vector>


namespace SM {
/*! \brief Virtual Hermitian operator matrix in p-p channel.
PPMatrix  represents an m-scheme matrix form of the operator \f$ \hat{H}=\sum_q \hat{T}_q \f$. The set of operators \f$ \hat{T}_q \quad q=0,1,2 \dots \f$ is stored as a vector of pointers to \ref sec_cioperator which should be incerted during initialization with <var>push_back</var> operator.

- The matrix is VIRTUAL, i.e. elements are computed on demand and it can be effeciently addressed with matrix vector multiplications. 

- Warning: change in operator results in change of the Matrix

  - Operators \f$ \hat{T}_q \f$ can be of different m-body structure. A Hamiltonian which contains one and two-body terms is a typical case of this matrix. 

We provide standard operations:
- Save to file with <var>Write</var> command SM::PPMatrix::Write();
- Multiplications with different <var>TimesVector</var> operations

\note non-Forward ordering will be ignored, locate assumes forward search and 
hermitian operator, SM::PPMatrix::TimesVector() also assumes hermicity 
\todo Combine Hermitian/Non-Hermitian
	*/
class PPMatrix: public std::vector<MBOperator*> {    
  public:
  uint_nbasis dim; ///< dimension of the matrix, identified as final states 
   int type; //0-normal 1-symmetric
   int SetSymmetric(){type=1;return 0;};
  int SetFull(){type=0; return 0;};
  const Many_Body_States &sti; ///<initial many-body states
  const Many_Body_States &stf; ///<final many-body states 
  /// Constructor, where initial and final set of states is different
  PPMatrix(Many_Body_States &_stf, Many_Body_States &_sti) : sti(_sti), stf(_stf) 
  {dim = stf.n;type=0;};
  /// Constructor for same initial and final states
  PPMatrix(Many_Body_States &_st) : sti(_st), stf(_st) 
  {dim = stf.n;type=1;}; 
   int TimesVector(double *yl, double *y, double *x, double shift);
   /*! \brief Matrix-vector multiplication, Hermitian (OpenMP)
   */
   int TimesVectorSymmetric(double *y, double *x, double shift); 
    int TimesVectorFull(double *y, double *x, double shift);
	int TimesVector(double *y, double* x, double shift){
   if (type==1) return TimesVectorSymmetric(y,x,shift); else return TimesVectorFull(y,x,shift);
  };
   /*! \brief get diatonal Hermitian (OpenMP)
   */
    int ReadDiagonal(double *y); 
	/*! \brief Separate diagonals for each of the constituent operators
	*/
	 int ReadDiagonals(double *y1, double *y2); 
	//<return  diagonal element defined by the state
	double DiagonalElement (int no, spsint *a) 
	{
	double retval=0.0;
	for (int io=0;io<(*this).size();io++) retval+=SM::ActDiagonalOperator((*(*this)[io]),no,a,1.0);
	return retval;
	};
    /*! \brief Deallocate componets
    This operator will clean all components in a vector and it assumes they were created with new. 
    It is a dangerous operator.
    */
    int DeAllocate (void)
	{
	for (int io=0;io<(*this).size();io++) {
    //((*this)[io]!=NULL)
    delete (*this)[io];
    }
	return 0;
	}; 
   /*! \brief Vector multiplication by the \f$ (TT^\dagger-S)|f\rangle \f$ operator, non-Hermitian (OpenMP). Fast version: via intermediate state. Version with forward operator
   is not possible. 
   both vectors are in ket-domain. 
   \test PPPower_test.cpp
   */
   int TimesVectorSqrNH(double *y, double *x, double shift); 
   /*! \brief Matrix Vector multiplication, \f$ (T-S)|i\rangle \f$ non-Hermitian (OpenMP)
   \test PPPower_test.cpp
   */
   int TimesVectorNH(double *y, double *x, double shift); 
   /*! \brief Non Hermitian k-th power Vector multiplication \f$ (T^k-S)|i\rangle \f$ (OpenMP)
   \test PPPower_test.cpp
   \note Similar approach can be used for power 2k
    */
   int TimesVectorPowNH(int k,double *y, double *x, double shift); 
 //  int TimesVector(double *y, double *x);
   /*! \brief Compute a square \f$ \langle i|T^\dagger T|i\rangle \f$ non-Hermitian operator
   \test PPPower_test.cpp
   */
   double Sqr(double *x);
   int Write(const char *filename); ///< Save operator to file [OpenMP] Hermitian
       //int ReadDiagonal(double *x);
   /*! Non Hermitian overlap non-Hermitian [OpenMP]
   \test SHLASF.cpp (computation of Spectroscopic Factors.)
    */
   double OverlapNH(const double *y, const double *x); 
    
  };
 int PPMatrix::TimesVectorSymmetric(double *y, double *x, double shift=0.0){
   if (sti.n==stf.n) for (uint_nbasis i=0;i<stf.n;i++) y[i]=-shift*x[i]; 
   else if (fabs(shift)<1E-8)  for (uint_nbasis i=0;i<stf.n;i++) y[i]=0.0;
   else FatalError("PPMatrix.TimesVector: Can not use shift if basis are not equal");
#pragma omp parallel default(shared)
 {
   uint_nbasis ket,locket;
   MBOperator Y(stf.N); //temporary operator
#pragma omp for schedule(dynamic) 
 for (uint_nbasis bra=0;bra<sti.n;bra++) {
 Y.clear();
 for (unsigned io=0;io<(*this).size();io++)
 SM::ActFermiOperator(Y,(*(*this)[io]),sti.N,sti.z[bra],1.0);
 locket=bra; //protect from out of range
 for (MBOperator::const_iterator iit=Y.begin();iit!=Y.end();iit++) {
   ket=locate(stf,iit->first,locket); //ket must be bigger than prevous ket
   if (ket==stf.n+1) continue;
   //Now I have matrix element iit-> second = H(bra,ket);
/* NB both operations must be done with atomic because other cpu writes y[ket] which may coniside with
y[bra] for current cpu.
*/
#pragma omp atomic
y[bra]+= (iit->second)*x[ket];
//critical causes some racing in rw   
 if (bra!=ket) 
{
#pragma omp atomic
 y[ket]+=(iit->second)*x[bra]; 
}
 locket=ket;
 }
 //cout<<bra<<" "<<nz<<endl;
 }
 } //end parallel region 
 return 0;
}

//Note this is an inverse operator multiplication

 int PPMatrix::TimesVectorFull(double *y, double *x, double shift=0.0){
   if (sti.n==stf.n) for (uint_nbasis i=0;i<stf.n;i++) y[i]=-shift*x[i]; 
   else if (fabs(shift)<1E-8)  for (uint_nbasis i=0;i<stf.n;i++) y[i]=0.0;
   else FatalError("PPMatrix.TimesVector: Can not use shift if basis are not equal");
  #ifdef CSTATUS
  av::status TM_status(sti.n,"TMF");
  #endif 

#pragma omp parallel default(shared)
 {
   uint_nbasis  ket;
   MBOperator Y(stf.N); //temporary operator
#pragma omp for schedule(dynamic) 
 for (uint_nbasis bra=0;bra<sti.n;bra++) {
 #ifdef CSTATUS
TM_status.ShowProgress(bra);
#endif 
 Y.clear();
 for (unsigned io=0;io<(*this).size();io++)
 SM::ActFermiOperator(Y,(*(*this)[io]),sti.N,sti.z[bra],1.0);
 
 uint_nbasis loket=0;
 for (MBOperator::iterator iit=Y.begin();iit!=Y.end();iit++) {
   ket=locate(stf,iit->first,loket); //ket may be anything
   if (ket==stf.n+1) continue;
   //Now I have matrix element iit-> second = H(bra,ket);   
 y[bra]+=(iit->second)*x[ket];
 loket=ket; //the search is done via loket because iit is sorted 
 }
 //cout<<bra<<" "<<nz<<endl;
 }
 } //end parallel region 
 return 0;
}

int PPMatrix::ReadDiagonal(double *y){
   if (sti.n!=stf.n) FatalError("PPMatrix.ReadDiagonal: Matrix must be square");
    #ifdef CSTATUS
   av::status RD_status(sti.n,"RD");
   #endif
#pragma omp parallel default(shared)
 {

#pragma omp for schedule(dynamic) 
 for (uint_nbasis bra=0;bra<sti.n;bra++) {
  #ifdef CSTATUS
 RD_status(bra);
 #endif
 y[bra]=DiagonalElement(sti.N,sti.z[bra]);
 //for (int io=0;io<(*this).size();io++) y[bra]+=SM::ActDiagonalOperator((*(*this)[io]),sti.N,sti.z[bra],1.0);
// SM::ActFermiOperator(Y,(*(*this)[io]),sti.N,sti.z[bra],1.0);
 //cout<<bra<<" "<<nz<<endl;
 }
 } //end parallel region 
 return 0;
}
//diagonals of each operator in matrix separately, y1 first, y2 second
int PPMatrix::ReadDiagonals(double *y1, double *y2){
   if (sti.n!=stf.n) FatalError("PPMatrix.ReadDiagonal: Matrix must be square");
    #ifdef CSTATUS
   av::status RD_status(sti.n,"RD");
   #endif
#pragma omp parallel default(shared)
 {

#pragma omp for schedule(dynamic) 
 for (uint_nbasis bra=0;bra<sti.n;bra++) {
 #ifdef CSTATUS
 RD_status(bra);
 #endif
 //y[bra]=DiagonalElement(sti.N,sti.z[bra]);
 y1[bra]=SM::ActDiagonalOperator((*(*this)[0]),sti.N,sti.z[bra],1.0);
 y2[bra]=SM::ActDiagonalOperator((*(*this)[1]),sti.N,sti.z[bra],1.0);
// SM::ActFermiOperator(Y,(*(*this)[io]),sti.N,sti.z[bra],1.0);
 //cout<<bra<<" "<<nz<<endl;
 }
 } //end parallel region 
 return 0;
}


/*! Non Hermitian version
need to use if operator is not hermitian: number of creation and annihilation not equal
\todo discuss "TimesVectorFull"
*/
int PPMatrix::TimesVectorNH(double *y, double *x, double shift=0.0){
  if (sti.n==stf.n) for (uint_nbasis i=0;i<stf.n;i++) y[i]=-shift*x[i]; 
  else if (fabs(shift)<1E-8)  for (uint_nbasis i=0;i<stf.n;i++) y[i]=0.0;
  else FatalError("PPMatrix.TimesVectorNH: Can not use shift if basis are not equal");
 // for (int i=0;i<stf.n;i++) y[i]=-shift*x[i]; //zero vector y   
#pragma omp parallel default(shared)
  { 
    uint_nbasis ket,locket;
    MBOperator Y(stf.N); //temporary operator 
#pragma omp for schedule(dynamic) 
  for (uint_nbasis bra=0;bra<sti.n;bra++) {
     Y.clear();
    for (unsigned io=0;io<(*this).size();io++)
      SM::ActFermiOperator(Y,(*(*this)[io]),sti.N,sti.z[bra],1.0);
    locket=0; //protect from out of range
    for (MBOperator::iterator iit=Y.begin();iit!=Y.end();iit++) {
      ket=locate(stf,iit->first,locket); //ket must be bigger than prevous ket
      if (ket==stf.n+1) continue;
  //if (locket >= ket) cerr<<"trouble :"<<locket<<" "<<ket<<endl;
   //Now I have matrix element iit-> second = H(bra,ket);
//This is important, to avoid racing condition we add y[bra]
#pragma omp atomic 
      y[ket]+=(iit->second)*x[bra]; 
      //y[bra]+=(iit->second)*x[ket]; //if basis the same we can do this
      locket=ket;
    }
 //cout<<bra<<" "<<nz<<endl;
  }
  }//end parallel region
  return 1;
}



double PPMatrix::OverlapNH(const double *y, const double *x){ 
  double htmp=0.0;
#pragma omp parallel default(shared) reduction(+:htmp)
  { 
    uint_nbasis ket,locket;
    MBOperator Y(stf.N); //temporary operator 
#pragma omp for schedule(dynamic) 
    for (uint_nbasis bra=0;bra<sti.n;bra++) {
      Y.clear();
      for (unsigned io=0;io<(*this).size();io++)
        SM::ActFermiOperator(Y,(*(*this)[io]),sti.N,sti.z[bra],1.0);
      locket=0; //protect from out of range
      for (MBOperator::iterator iit=Y.begin();iit!=Y.end();iit++) {
        ket=locate(stf,iit->first,locket); //ket must be bigger than prevous ket
        if (ket==stf.n+1) continue;
  //if (locket >= ket) cerr<<"trouble :"<<locket<<" "<<ket<<endl;
   //Now I have matrix element iit-> second = H(bra,ket);
//This is important, to avoid racing condition we add y[bra]
//#pragma omp atomic does not need pragma
        htmp+=y[ket]*(iit->second)*x[bra];
      //y[bra]+=(iit->second)*x[ket]; //if basis the same we can do this
        locket=ket;
      }
 //cout<<bra<<" "<<nz<<endl;
    }
  }//end parallel region
  return htmp;
}



int PPMatrix::TimesVectorPowNH(int k, double *y, double *x, double shift=0.0){
  //for (int i=0;i<stf.n;i++) y[i]=-shift*x[i]; //zero vector y   
  if (sti.n==stf.n) for (uint_nbasis i=0;i<stf.n;i++) y[i]=-shift*x[i]; 
  else if (fabs(shift)<1E-8)  for (int i=0;i<stf.n;i++) y[i]=0.0;
  else FatalError("PPMatrix.TimesVectorPowNH: Can not use shift if basis are not equal");
#pragma omp parallel default(shared)
  { 
    uint_nbasis ket,locket;
    MBOperator Z(stf.N), ZP(stf.N);
    MBOperator *Y=&Z; //temporary operator 
    MBOperator *YP=&ZP; //temporary operator 
#pragma omp for schedule(dynamic) 
    for (uint_nbasis bra=0;bra<sti.n;bra++) {
      Y->clear();
      (*Y)[const_cast<const spsint*>(sti.z[bra])]=1.0; //first Y
      /*here we compute power of the operator:
      we enter with original y and iterate it k times 
      */
      for (int p=0;p<k;p++) {
      YP->clear();
      for (MBOperator::iterator iit=Y->begin();iit!=Y->end();iit++) 
      for (unsigned io=0;io<(*this).size();io++)
        SM::ActFermiOperator((*YP),(*(*this)[io]),sti.N,iit->first,iit->second);
      
      //swap pointers;
      MBOperator *Ytmp=Y; Y=YP; YP=Ytmp;
      //Y.clear();
      //Y+=YP; //new Y
      }
      
      
      locket=0; //protect from out of range
      for (MBOperator::iterator iit=Y->begin();iit!=Y->end();iit++) {
        ket=locate(stf,iit->first,locket); //ket must be bigger than prevous ket
        if (ket==stf.n+1) continue;
  //if (locket >= ket) cerr<<"trouble :"<<locket<<" "<<ket<<endl;
   //Now I have matrix element iit-> second = H(bra,ket);
//This is important, to avoid racing condition we add y[bra]
#pragma omp atomic 
        y[ket]+=(iit->second)*x[bra]; //we write transposed matrix
        locket=ket;
      }
 //cout<<bra<<" "<<nz<<endl;
    }
  }//end parallel region
  return 1;
}



int PPMatrix::TimesVectorSqrNH(double *y, double *x, double shift=0.0){
 // for (int i=0;i<stf.n;i++) y[i]=-shift*x[i]; //zero vector y   
  if (sti.n==stf.n) for (uint_nbasis i=0;i<stf.n;i++) y[i]=-shift*x[i]; 
  else if (fabs(shift)<1E-8)  for (uint_nbasis i=0;i<stf.n;i++) y[i]=0.0;
  else FatalError("PPMatrix.TimesVectorSqrNH: Can not use shift if basis are not equal");
#pragma omp parallel default(shared)
  { 
    uint_nbasis ket,locket;
    MBOperator Y(stf.N); //temporary operator 
#pragma omp for schedule(dynamic) 
    for (uint_nbasis bra=0;bra<sti.n;bra++) {
      Y.clear();
      for (unsigned io=0;io<(*this).size();io++)
        SM::ActFermiOperator(Y,(*(*this)[io]),sti.N,sti.z[bra],1.0);
      locket=0; //protect from out of range
      unsigned nz=0;
      uint_nbasis *ind= new uint_nbasis[Y.size()];
      double *elem=new double [Y.size()];
      //temporarily store data
      for (MBOperator::iterator iit=Y.begin();iit!=Y.end();iit++) {
        ket=locate(stf,iit->first,locket); //ket must be bigger than prevous ket
        if (ket==stf.n+1) continue;
        ind[nz]=ket;
        elem[nz]=iit->second;
        nz++;
        locket=ket;
      }
      //do the actual computation
     // cout<<bra<<" "<<nz<<endl;
      //Pause();

      for (unsigned ki=0;ki<nz;ki++) 
        for (unsigned kf=0;kf<nz;kf++) 
#pragma omp atomic  
          y[ind[kf]]+=elem[ki]*elem[kf]*x[ind[ki]];
 //cout<<bra<<" "<<nz<<endl;
      delete [] ind;
      delete [] elem;
    }
  }//end parallel
  return 1;
}





double PPMatrix::Sqr(double *x){
    MBOperator Y(stf.N); //temporary operator 
    for (uint_nbasis bra=0;bra<sti.n;bra++) {
      for (uint_nbasis io=0;io<(*this).size();io++)
        SM::ActFermiOperator(Y,(*(*this)[io]),sti.N,sti.z[bra],x[bra]);
      //Open MP need reduction
    }
  return SM::Sqr(Y);
}

 int PPMatrix::Write(const char *filename){
   ofstream outf;
   outf.open(filename, ios::binary); //start my opening file
   outf.write((char*)(&sti.n),sizeof(uint_nbasis)); //we have saved n
    av::status WH_status(sti.n,"WH");
#pragma omp parallel default(shared)
   { 
 unsigned nz; //number of non-zero elements;
 double *elem = new double[stf.n];
 uint_nbasis *ind=new uint_nbasis[stf.n];
 uint_nbasis ket,locket;
 MBOperator Y(stf.N); //temporary operator
#pragma omp for ordered schedule(dynamic) 
 for (uint_nbasis bra=0;bra<sti.n;bra++) {
 //StatusPercent(bra, sti.n);
  WH_status(bra);
 Y.clear();
 //io is the component in pp matrix.
 for (unsigned io=0;io<(*this).size();io++) {
 //if(sti.N>=(((*this)[io])->rank-Y.rank)) //the check is automatically done in ActFermi
 SM::ActFermiOperator(Y,(*(*this)[io]),sti.N,sti.z[bra],1.0);
 }
 nz=0;
 locket=bra; //protect from out of range
 for (MBOperator::iterator iit=Y.begin();iit!=Y.end();iit++) {
   if (type) ket=locate(stf,iit->first,locket); //ket must be bigger than prevous ket
   else ket=locate(stf,iit->first);
   if (ket==stf.n+1) continue;
   ind[nz]=locket=ket;
  // cout<<bra<<" "<<ket<<" "<<(iit->second)<<endl;
   elem[nz]=(iit->second);
   nz++; 
 }
 //cout<<bra<<" "<<nz<<endl;

#pragma omp ordered 
 {
 // cerr<<bra<<" "<<nz<<" ";
 // for (int i=0;i<nz;i++) cout<<ind[i]<<" ";
 // cout<<endl;
 outf.write((char*)(&(nz)),sizeof(unsigned)); //we saved nz
 outf.write((char*)(ind),nz*sizeof(uint_nbasis)); //saved array of ind
 outf.write((char*)(elem),nz*sizeof(double)); //saved array of elem
 }
 
 //cout<<Y<<endl;
 }
 //cout<<stf.N<<endl;
 delete [] elem;
 delete [] ind;
   } //end parallel region
 return 1;
}



#ifdef __TAV_SPARSE_CXX__
    //save into sparse matrix
    int PPMatrix2SparseMatrix(tav::SparseMatrix<double> &X, PPMatrix &HH)
    {
        X.type=1; //for symmetric
#ifdef CSTATUS
        av::status VH_status(HH.sti.n,"VH");
#endif
#pragma omp parallel default(shared)
        {
            unsigned nz; //number of non-zero elements;
            double *elem = new double[HH.stf.n];
            uint_nbasis *ind=new uint_nbasis[HH.stf.n];
            uint_nbasis ket,locket;
            MBOperator Y(HH.stf.N); //temporary operator
#pragma omp for schedule(dynamic)
            for (uint_nbasis bra=0;bra<HH.sti.n;bra++)
            {
#ifdef CSTATUS
                VH_status(bra);
#endif
                //cerr<<bra<<" "<<HH.sti.n<<endl;
                Y.clear();
                for (unsigned io=0;io<HH.size();io++)
                    SM::ActFermiOperator(Y,(*HH[io]),HH.sti.N,HH.sti.z[bra],1.0);
                nz=0;
                locket=bra; //protect from out of range
                for (MBOperator::iterator iit=Y.begin();iit!=Y.end();iit++)
                {
                    ket=locate(HH.stf,iit->first,locket); //ket must be bigger than prevous ket
                    if (ket==HH.stf.n+1) continue;
                    ind[nz]=locket=ket;
                    elem[nz]=(iit->second);
                    nz++;
                }
                //cout<<bra<<" "<<nz<<endl;
                
                //#pragma omp critical //note potential racing condition in Sync operator
                X.WriteLine(bra, nz, ind, elem);
            }
            //cout<<stf.N<<endl;
            delete [] elem;
            delete [] ind;
        } //end parallel region
        return 1;
    }
#endif 
    
#ifdef  __Tav__CSBMatrix__
    
    //save into sparse matrix
    template<class ubig_t = uint32_t, class usmall_t = uint16_t>
    int PPMatrix2CSBRows(tav::SparseMatrixCSB<double,ubig_t,usmall_t> &X, PPMatrix &HH)
    {
        X.isSymmetric=true; //for symmetric
#ifdef CSTATUS
        av::status VH_status(HH.sti.n,"VH");
#endif
#pragma omp parallel default(shared)
        {
            ubig_t nz; //number of non-zero elements;
            double *elem = new double[HH.stf.n];
            ubig_t *ind=new ubig_t[HH.stf.n];
            uint_nbasis ket,locket;
            MBOperator Y(HH.stf.N); //temporary operator
#pragma omp for schedule(static) ordered
            for (uint_nbasis bra=0;bra<HH.sti.n;bra++)
            {
#ifdef CSTATUS
                VH_status(bra);
#endif
                //cerr<<bra<<" "<<HH.sti.n<<endl;
                Y.clear();
                for (unsigned io=0;io<HH.size();io++)
                    SM::ActFermiOperator(Y,(*HH[io]),HH.sti.N,HH.sti.z[bra],1.0);
                nz=0;
                locket=bra; //protect from out of range
                for (MBOperator::iterator iit=Y.begin();iit!=Y.end();iit++)
                {
                    ket=locate(HH.stf,iit->first,locket); //ket must be bigger than prevous ket
                    if (ket==HH.stf.n+1) continue;
                    locket=ket;
                    ind[nz]= ubig_t(ket);
                    elem[nz]=(iit->second);
                    nz++;
                }
                //cout<<bra<<" "<<nz<<endl;
                
                //#pragma omp critical //note potential racing condition in Sync operator
                X.WriteRowMove(bra, elem, ind, nz);
            }
            //cout<<stf.N<<endl;
            delete [] elem;
            delete [] ind;
        } //end parallel region
        
        X.MakeDiagonals();
        
        return 1;
    }
    
    template<class ubig_t = uint32_t, class usmall_t = uint16_t>
    int PPMatrix2CSBBlk(tav::SparseMatrixCSB<double,ubig_t,usmall_t> &X, PPMatrix &HH)
    {
        X.isSymmetric=true; //for symmetric#endif

        //assume everything is correctly sized.
        //re-save for ease.
        usmall_t BlockRow = X.BlockRow, BlockCol = X.BlockCol;
        usmall_t numBlockRow = X.numBlockRow, numBlockCol = X.numBlockCol;
        
    #ifdef CSTATUS
            av::status VH_status(numBlockRow,"VH");
    #endif
#pragma omp parallel for
        for (usmall_t row = 0; row < numBlockRow; ++row)
        {
            std::vector<double>* elem = new std::vector<double> [numBlockCol];
            std::vector<usmall_t>* indCol = new std::vector<usmall_t> [numBlockCol];
            std::vector<usmall_t>* indRow = new std::vector<usmall_t> [numBlockCol];
            for (usmall_t i = 0; i < numBlockCol; ++i)
            {
                X.Blocks[i][i].isDiagonal = true;
                //reserving makes it faster (duh), the factor of 10 
                //is arbitrary. but gave better results
                //than some other values I don't remember right now
                //feel free to experiment :) 
                elem[i].reserve(BlockRow * BlockCol / 1000);
                indCol[i].reserve(BlockRow * BlockCol / 1000);
                indRow[i].reserve(BlockRow * BlockCol / 1000);
            }
            
            ubig_t* nnz = new ubig_t[numBlockCol];
            memset(nnz,0,numBlockCol*sizeof(ubig_t));//zero all non-zero matrix element counters
            ubig_t finalbra = ( (row + 1) * BlockRow < HH.sti.n ? (row+1)*BlockRow : HH.sti.n);
            
            for (ubig_t bra=row * BlockRow; bra < finalbra; bra++)
            {
                ubig_t locket,ket;
                usmall_t blkc, blkr;

                
                MBOperator Y(HH.stf.N);
                for (unsigned io=0;io<HH.size();io++)
                    SM::ActFermiOperator(Y,(*HH[io]),HH.sti.N,HH.sti.z[bra],1.0);
                //StripZeros(Y);
                
                locket=bra; //protect from out of range
                for (MBOperator::iterator iit=Y.begin();iit!=Y.end();iit++)
                {
                    ket=locate(HH.stf,iit->first,locket); //ket must be bigger than prevous ket
                    if (ket==HH.stf.n+1) continue;
                    if (std::abs(iit->second) < 1.e-14) continue;
                    locket=ket;
                    blkc = ket / BlockCol;
                    indCol[blkc].push_back( ket % BlockCol );
                    indRow[blkc].push_back( bra % BlockRow );
                    elem[blkc].push_back( iit->second );
                    nnz[blkc]++;
                }
            }
            for (usmall_t j = 0; j < numBlockCol; ++j)
            {
                if (nnz[j] > 0)
                {
                    X.Blocks[row][j].WriteFullBlockMove(&(indRow[j][0]),&(indCol[j][0]),&(elem[j][0]),nnz[j]);
                }
            }
    #ifdef CSTATUS
                VH_status(row);
    #endif
            
            delete[] nnz;
//            for (usmall_t j = 0; j < numBlockCol; ++j)
//            {
//                delete[] elem[j];
//                delete[] indCol[j];
//                delete[] indRow[j];
//            }
            delete[] elem;
            delete[] indCol;
            delete[] indRow;
            
        }//end rows
        
        X.MakeDiagonals();
        return 1;
    }
    
#endif
#ifdef EIGEN_WORLD_VERSION
    int PPMatrix2MatrixXd(Eigen::MatrixXd &X, PPMatrix &HH)
    {
        if (HH.sti.N == HH.stf.N && HH.stf.N == 0)
        {
            X=Eigen::MatrixXd::Zero(1,1);//There can only be 1 state in this casse.
            
            std::cerr << "MBS is for closed core (no particles) returning 0 matrix." << std::endl;
            X(0,0) = 0.;
            return -1;
        }
        X=Eigen::MatrixXd::Zero(HH.stf.n,HH.sti.n);
#ifdef CSTATUS
        av::status VH_status(HH.sti.n,"VH");
#endif
#pragma omp parallel default(shared)
        {
            unsigned nz; //number of non-zero elements;
//            double *elem = new double[HH.stf.n];
//            uint_nbasis *ind=new uint_nbasis[HH.stf.n];
            uint_nbasis ket,locket;
            MBOperator Y(HH.stf.N); //temporary operator
#pragma omp for schedule(dynamic)
            for (uint_nbasis bra=0;bra<HH.sti.n;bra++)
            {
#ifdef CSTATUS
                VH_status(bra);
#endif
                //cerr<<bra<<" "<<HH.sti.n<<endl;
                Y.clear();
                for (unsigned io=0;io<HH.size();io++)
                    SM::ActFermiOperator(Y,(*HH[io]),HH.sti.N,HH.sti.z[bra],1.0);
                nz=0;
                locket=bra; //protect from out of range
                for (MBOperator::iterator iit=Y.begin();iit!=Y.end();iit++) {
                    ket=locate(HH.stf,iit->first,locket); //ket must be bigger than prevous ket
                    if (ket==HH.stf.n+1) continue;
//                    ind[nz]=locket=ket;
                    locket=ket;
                    X(locket,bra)=X(bra,locket)=iit->second;
//                    elem[nz]=(iit->second);
                    nz++;
                }
                //cout<<bra<<" "<<nz<<endl;
                
                //#pragma omp critical //note potential racing condition in Sync operator
//                X.WriteLine(bra, nz, ind, elem);
            }
            //cout<<stf.N<<endl;
//            delete [] elem;
//            delete [] ind;
        } //end parallel region
        return 1;
    }
#endif
/*
int SHLMatrix::TimesVector(double *y, double *x){
  return TimesVector(y,x,0.0);
}

int SHLMatrix::TimesVector(double *yl, double *y, double *x, double shift){
double *elem=new double [st.n];  //value of non-zero element
int *ind= new int [st.n]; //location of non-zero element
int nz; //how many non zero elements
INF(cerr<<"Generating Hamiltonian Matrix     "<<endl;)


 //actual multiplication procedure
 for (int i=0;i<st.n;i++) y[i]=0.0; //zero vector y
 for (int i=0;i<st.n;i++) {
   nz= ActFermiOperatorX(elem, ind, i, st, ESP, VSP );
  
   //to avoid new and delete we keep space for full dimension n
   //NOTE WE COULD REDUCE MEMORY USE BY SIZING ARRAY TO NZ
   if (nz==0) y[i]=-shift*x[i]; // if there are no elements
   else if (ind[0]==i) elem[0]-=shift; //shift diagonal
                for (int k=0;k<nz;k++) {
		  y[i]+=elem[k]*x[ind[k]];
                  if(i!=ind[k]) y[ind[k]]+=elem[k]*x[i];}
 
    
 }
 for (int i=0;i<st.n;i++) yl[i]=y[i];
 delete [] elem;
 delete [] ind;
 return 1;
}

*/
} //SM space name
#endif

