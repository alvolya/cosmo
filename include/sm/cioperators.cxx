/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file cioperators.cxx
\brief Configuration Interaction Many-Body Operators

*/
/*
06/10/2009 New operator SPSqrM produces non-symmetrized version
10/25/2007 created as a more advanced version of operators
 */
 
 
 
 
#ifndef __CIOPERATORS_CXX__
#define __CIOPERATORS_CXX__
#include <cmath> //needed for sqrt
#include "operators.cxx"
#include <tav/stensor.cxx>
#include <vector>
using std::vector;
#define _ZERO_TOLERANCE_ 1E-8
using tav::STensor;
/// many-body operator coded as sparse tensor 
typedef tav::STensor<spsint, double> MBOperator; 
namespace SM {
/*! \brief Operation \f$ |f\rangle+=\hat{Q}\, c|i\rangle \f$
*/
int ActFermiOperator(MBOperator &fstate, ///< Final state (output), rank=number of particles 
		 const MBOperator &QQ,   ///< n-body operator, rank=na+(rank-na)
		 int NN,	///< number of particles in initial state
		 spsint *istate,     ///< pointer to initial state
		 double C=1.0        ///< amplitude of the initial state, default [1]
		) 
{
	int na=((QQ.rank+NN-fstate.rank)>>1); //number of annihilation operators
	//On=na+nc; Nf=Ni-na+nc
    if (NN<na) return 1; //There are more annihilaiton operators then particles
	int qn; //temporary integer for particle number in intermediate state
	spsint *newspstate=new spsint [fstate.rank]; //final state
	spsint *newspstate1=new spsint [NN-na]; //dimension of intermediate state, after annihilaiton
	bool addphase=false;
	//MBOperator::iterator iit; //iterator to scan through oper
	for (auto iit=QQ.begin(); iit!=QQ.end(); iit++) {
		spsint *oper=(*iit).first; //pointer to operator
		addphase=false; //initialize phase
	if(!Annihilation(addphase,qn,newspstate1,na,oper,NN,istate)) continue;
	if(!Creation(addphase,qn,newspstate,QQ.rank-na,(oper+na),qn,newspstate1)) continue;
        ///\todo The incertion can be accelerated given that states are sorted
	fstate[newspstate]+=(addphase? -(*iit).second * C: (*iit).second * C); 
	}
	delete [] newspstate; 
	delete [] newspstate1;
	return 0;
}

/*! \brief Operation \f$ |f\rangle+=\hat{Q}\, c|i\rangle \f$
We only select diagonal, so final state is the same
*/
double ActDiagonalOperator(
		 MBOperator &QQ,   ///< n-body operator, rank=na+(rank-na) 
		 int NN,	///< number of particles in initial state
		 spsint *istate,     ///< pointer to initial state
		 double C=1.0        ///< amplitude of the initial state, default [1]
		) 
{
	double retvalue=0.0;
	int na=((QQ.rank)>>1); //number of creation operators
	//On=na+nc; Nf=Ni-na+nc
	int qn; //temporary integer for particle number in intermediate state
	spsint *newspstate=new spsint [NN]; //temporary operators
	bool addphase=false;
	MBOperator::iterator iit; //iterator to scan through oper
	for (iit=QQ.begin(); iit!=QQ.end(); iit++) {	    
		spsint *oper=(*iit).first; //pointer to operator
		if (CompareGeneral(na, oper, oper+na)!=0) continue; //if operators are not diagonal continue
		addphase=false; //initialize phase
	if(!Annihilation(addphase,qn,newspstate,na,oper,NN,istate)) continue;
	//if(!Creation(addphase,qn,newspstate,QQ.rank-na,(oper+na),qn,newspstate1)) continue;
        ///\todo The incertion can be accelerated given that states are sorted
	retvalue+=(*iit).second * C; 
	}
	delete [] newspstate; 
	return retvalue;
}


/*! \brief Operation \f$ \hat{F}+=\hat{Q}\, \hat{I} \f$ All operators are SM::MBOperator
 */
int Product(MBOperator &fstate, ///< Final state (output), rank=number of particles 
                     const MBOperator &QQ,   ///< n-body operator, rank=na+(rank-na)
                     const MBOperator &istate, ///< Initial state (input), rank=number of particles
                    double C = 1.0)
{
  for (auto iit=istate.begin();iit!=istate.end();iit++) //scan over initial
    SM::ActFermiOperator(fstate, QQ, istate.rank,(iit->first),iit->second * C);
  
  return 1;
}

    /*! \brief Parallel version of SM::Product */
    int ParallelProduct(MBOperator& fstate,
                        const MBOperator& QQ,
                        const MBOperator& istate,
                        double C=1.0)
    {
        std::vector<MBOperator::const_iterator> its(istate.size(),istate.end());
        size_t ixt = 0;
        for (auto it = istate.begin(); it != istate.end(); ++it,++ixt)
                its[ixt] = it;
        //std::vector<MBOperator> results(istate.size(),MBOperator(fstate.rank));
#pragma omp parallel for
        for (int i = 0; i < its.size(); ++i)
        {
        MBOperator x(fstate.rank);
 //       std::cout << "OP: " << i << std::endl;     
           SM::ActFermiOperator(x, QQ, istate.rank,its[i]->first, its[i]->second * C);
#pragma omp critical
            {
                fstate += x;
            }
        }
        
//        for (size_t i = 0; i < results.size(); ++i)
//            fstate += results[i];
        return 1;
    }

//this is already present in STensor
double Sqr(const MBOperator &QQ) {return tav::Sqr(QQ);}
///*! \brief Find Square of the State \f$ \langle 0|T T^\dagger|0\rangle \f$
// */
//double Sqr(const MBOperator &QQ   ///< n-body operator, rank=na+(rank-na)
//             )
//{
//  double htmp=0.0;
//  for (auto iit=QQ.begin();iit!=QQ.end();iit++) //scan over initial
//    htmp+=(iit->second)*(iit->second);
//  return htmp;
//}

/*! \brief Find Normalization of the State
 */
double Abs(const MBOperator &QQ   ///< n-body operator, rank=na+(rank-na)
          ) 
{
   return sqrt(tav::Sqr(QQ));
}

/*! \brief Dot product of two states
 */
double Dot(const MBOperator &bra,    ///< n-body operator,
const MBOperator &ket //ket operator
)
{
double x=0.0;
if (bra.rank!=ket.rank) return 0.0;
for (auto iit=bra.begin();iit!=bra.end();iit++)
{
auto iik=ket.find((iit->first));
if (iik==ket.end()) continue;
else
x+=(iit->second)*(iik->second);

}
return x;
}

/*! \brief Operation \f$ \hat{F}=\hat{Q}^\dagger \hat{P} \f$
\todo discuss ordering and conjugation
*/
int ConjugateProduct(MBOperator &F, ///< Final state (output), rank=number of particles 
                     MBOperator &Q,   ///< read as annihilation ...
                     MBOperator &P ///< Initial state (input), rank=number of particles
                    ) 
{
  spsint *a=new spsint[F.rank]; //maybe a problem if F.rank!=P.rank+Q.rank
  for (MBOperator::iterator ip=P.begin();ip!=P.end();ip++) //scan over initial
    for (MBOperator::iterator iq=Q.begin();iq!=Q.end();iq++) {
    for (int i=0;i<P.rank;i++) a[i]=(ip->first)[i];
    for (int i=0;i<Q.rank;i++) a[i+P.rank]=(iq->first)[i];
    F[a]=+(ip->second)*(iq->second); //note in principal first must be conjugated
    }
  delete []a;
  return 1;
}


		   /*! \brief Square particle hole operator leading to 2-body and 1-body M^+ M
The operator is stored in one and two-body parts. If initial is not zero than we add. The output is not striped for 
hermicity
*/		 
 int SPSqrM(tav::STensor<spsint,double> &JSP,  
            tav::STensor<spsint,double> &VSP, 
            tav::STensor<spsint,double> &JO 
           ) {
//we add 
    //JSP=0.0;
    //VSP.clear();

             spsint aa[4];
             tav::STensor<spsint,double>::iterator iit;
             tav::STensor<spsint,double>::iterator jjt;
             for (iit=JO.begin();iit!=JO.end();iit++) {
     //oprator a+(i1) a(i0)
               for (jjt=JO.begin();jjt!=JO.end();jjt++){
			   double htmp=(*iit).second*(*jjt).second;
			   if (fabs(htmp)<_ZERO_TOLERANCE_) continue;
       /* total oprator [a+(j1) a(j0)]^+  a+(i1) a(i0)=
                 a^+(j0) a(j1) a^+(i1) a(i0)= delta(i1,j1) a(j0)^2 a(i0) - a^+(j0)  a^+(i1) a(j1) a(i0)=
                 ...delta(i1,j1) a(j0)^+ a(i0) -[a(i1)  a(j0)]^+ a(j1) a(i0)
       */

                 if ((*jjt).first[1]==(*iit).first[1])
                   JSP((*jjt).first[0],(*iit).first[0])+=htmp;
       //the format here is JSP[i][j]a^(i) a(j);
                  
                 aa[0]=(*jjt).first[1]; //j1
				 aa[1]=(*iit).first[0]; //i0
				 aa[2]=(*jjt).first[0]; //j0
				 aa[3]=(*iit).first[1]; //i1 //pay attention to extra minus in equation
				 
				 if (aa[0]>aa[1]) {tav::Swap(aa[0],aa[1]); htmp=-htmp;}
				 if (aa[2]>aa[3]) {tav::Swap(aa[2],aa[3]); htmp=-htmp;}
                   VSP[aa]+=htmp; 
                 }
        
               }
// delete [] aa;
            //SM::StripHermitianOperator(JSP);
			//SM::StripHermitianOperator(VSP);
               return 1;
           };

/*! \brief Converted a forward acting operator to a symmetric full operator by dublicating matrix elements
           */
int Hermitian2FullOperator(MBOperator &VX){
  int na=VX.rank>>1; ; //get number of operators
  spsint *a =new spsint [VX.rank]; //new array
  MBOperator::iterator iit; //iterator to scan through oper
  for (iit=VX.begin(); iit!=VX.end(); iit++) {
   
    spsint *x=iit->first; 
    spsint *y=x+na;
    for (int i=0;i<na;i++) { 
      a[i]=y[i]; a[i+na]=x[i];  
    }
    //check if this is diagonal
    bool flag=false;
    for (int i=0;i<na;i++) if (a[i]!=a[i+na]) flag=true;
    if (flag) VX[a]=iit->second; //make sure this is not diagonal
    /*
    for (int i=0;i<VX.rank;i++) cout<<int(x[i])<<" ";
    cout<<" : ";
    for (int i=0;i<VX.rank;i++) cout<<int(a[i])<<" ";
    cout<<" flag= "<<flag<<endl;
    Pause();
    }
    */
    /// \todo there should be VX[a]+= command but generation of forward operator from x is not correct  fix: SM::MakeVSPO(VSPX,sps,x);
} 
return 0;
}

/*! \brief Converted an arbitrary operator to a forward acting hermitian
 */
int StripHermitianOperator(MBOperator &VX){
  int na=VX.rank>>1; ; //get number of operators
  spsint *a =new spsint [VX.rank]; //new array
  MBOperator::iterator iit; //iterator to scan through oper
  tav::ArrayComparison<spsint> XYZ(na);
  for (iit=VX.begin(); iit!=VX.end(); iit++) {
   
  //  spsint *x=iit->first; 
   // spsint *y=x+na;
    //the following approach is safe but not fast
    if (XYZ((iit->first)+na,iit->first)) {VX.erase(iit); iit=VX.begin(); }//reset pointer
  } 
  return 1;
}

    int StripHermitianOperator(MBOperator &Temp, const MBOperator &VX)
    {
        //MBOperator Temp(VX.rank);
        Temp.clear();
        Temp.SetRank(VX.rank);
        int na=VX.rank>>1; ; //get number of operators
        spsint *a =new spsint [VX.rank]; //new array
        //MBOperator::iterator iit; //iterator to scan through oper
        tav::ArrayComparison<spsint> XYZ(na);
        for (auto iit=VX.begin(); iit!=VX.end(); iit++)
        {
            if (std::abs(iit->second)<1.e-12) continue;
            //  spsint *x=iit->first;
            // spsint *y=x+na;
            //the following approach is safe but not fast
            if (!(XYZ((iit->first)+na,iit->first))) {Temp[iit->first]=iit->second;}//reset pointer
            
        }
        //VX.clear();
       // VX=Temp;
        return 0;

        
    }
    /*! \brief The goal is to strip the operator of all zeros and non hermitian elements.
     it stores to new operator and then replaces it.
     */
     //forer name StripWithCopy
    int StripHermitianOperatorWithCopy(MBOperator &VX)
    {
        MBOperator Temp(VX.rank);
StripHermitianOperator(Temp, VX);
        VX.clear();
        VX=Temp;
        return 1;

        
    }
/*! \brief Converted an arbitrary operator to a forward acting hermitian
 */
int Conjugate(MBOperator &XV, MBOperator &VX, 
              int na ///<number of annihilation operators
             ){
  spsint *a =new spsint [VX.rank]; //new array
  MBOperator::iterator iit; //iterator to scan through oper
  for (iit=VX.begin(); iit!=VX.end(); iit++) {
     for (int i=0;i<VX.rank-na;i++) a[i]=(iit->first)[i+na];
     for (int i=VX.rank-na;i<VX.rank;i++) a[i]=(iit->first)[i-VX.rank+na];
     XV[a]+=iit->second;
  }
  delete [] a;
  return 1;
}


/*! \brief This function takes an MBOperator and converts the key from a first set of sps to another one.
 Result is saved in a diffrent operator since keys are const by nature and const_casting is not a good
 programming practice according to people I don't know on the internet.
 */
int ConvertMBOperatorKeys(MBOperator& result,MBOperator& op, vector<spsint> mappings, double tol=_ZERO_TOLERANCE_)
{
    spsint k;
    spsint *key = new spsint [op.rank];
    bool saveflag=true;
    
    for (MBOperator::iterator itx = op.begin(); itx != op.end(); itx++)
    {
    //cerr<<"op "<<(itx->second)<<endl;
        if (std::abs(itx->second) < tol) continue;//if value too small there is no point in saving it.
        saveflag = true;
        for (size_t i = 0; i < op.rank; i++)
        {
            k = itx->first[i];
            if (mappings[k] == spsint(-1))
            {
                saveflag=false;
                break;
            }
            key[i] = mappings[k];
        }
        if (saveflag)
        {
//            for (int i =0;i<op.rank;i++)
//                std::cout << key[i] << " ";
//            std::cout << itx->second << std::endl;
//            std::sort(key,key+op.rank);//no need usually but just in case. (comment out maybe?)
            int phase = 1;
            for (int i = 0; i < op.rank; ++i)
            {
                for (int j = i+1; j <op.rank; ++j)
                    if (key[i] > key[j])
                    {
                        phase = -phase;
                        std::swap(key[i],key[j]);
                    }
            }
            result[key] = itx->second * phase;
        }
    }
    delete[] key;
    return result.size();
}
    
    /*! \brief Find minimum and maximum key in operator
 */
int MinMaxKeys(spsint &min, spsint &max,
    const MBOperator &VX //input operator
             ){
    if (VX.size()==0) return 1; //return one if operator is empty
   spsint lst=VX.rank-1; //set location of last
   auto iit=VX.begin(); //do first step manually
  min=(iit->first)[0];
  max=(iit->first)[lst];
  iit++; //manually increment
  for (; iit!=VX.end(); iit++) {
   if (min>(iit->first)[0]) min=(iit->first)[0];
   if (max<(iit->first)[lst]) max=(iit->first)[lst];
  }
  return 0;
}


} //end namespace
#endif 

