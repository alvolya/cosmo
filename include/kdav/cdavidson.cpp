
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <ctype.h>
#include <cstring>
using std::strcpy;
#include "cdavidson.h"
typedef unsigned long long uint_t_ee_file; //unsigned interger used in EE file

/*This is an initialization function of davidson process targeting a sequence of eivenvalues
 */
Davidson::Davidson(
                   const MatB::MatrixBase *diag, ///< matrix,
                   ntyp eigvalmax, ///< Highest wanted eigenvalue 
                   ntyp eigvalmin) ///< lowest wanted eigenvalue
:ilow(eigvalmin),n(diag->Size()),ihigh(eigvalmax),tab(false),diag(diag),iiwsz(0),lim(0),irwsz(0),neig(0),ierr(0),niv(0),nmv(0),nume(0),numemax(0),hiend(false),Init(false),Init2(false),b(),ab(),twosvec(),smt(),eigs(),scratch1(),old(),icv(),scra2(),scra3(),critr(0),critc(0),crite(0),ortho(0),sum(0),maxiter(0),iselec(NULL),silence(false),mblock(0),limmax(0),norm(0),loop(0),inf(false)

{
    
//    Davidson::diag = diag; //copy pointer to internal variable
//    n = diag->Size(); //size of the matirx
//    tab=false;
//
//    ilow = eigvalmin;
//    ihigh = eigvalmax;
    if (!ilow) ilow=ihigh;

    ntyp Nume;
    if ((ilow+ihigh-1)>n) {
        Nume=n-ilow+1;
    } else Nume=ihigh; //set Nume to the number of the heighest eigenvalue, eigenvalues counted here as 1,2...eigvalmax.
    SetDefault(Nume); //setup 

    iselec = new ntyp[limmax]; //mona by poprawic na liczbe poszukiwanych wart. wasnych//can be improved on the number of wanted worth. important
}

Davidson::Davidson(const MatB::MatrixBase *diag, ntyp eigval[])
:ilow(0),n(diag->Size()),ihigh(0),tab(false),diag(diag),iiwsz(0),lim(0),irwsz(0),neig(0),ierr(0),niv(0),nmv(0),nume(0),numemax(0),hiend(false),Init(false),Init2(false),b(),ab(),twosvec(),smt(),eigs(),scratch1(),old(),icv(),scra2(),scra3(),critr(0),critc(0),crite(0),ortho(0),sum(0),maxiter(0),iselec(NULL),silence(false),mblock(0),limmax(0),norm(0),loop(0),inf(false)

{
//    Davidson::diag = diag;
//    n = diag->Size();
    tab=true;

    ilow = 0;
    ntyp eigvallow=eigval[0], eigvalhigh=eigval[0], Nume;
    for (ntyp i=1; eigval[i]>0; i++) {
    	eigvallow=MIN(eigvallow, eigval[i]);
        eigvalhigh=MAX(eigvalhigh, eigval[i]);
    }
    if ((eigvallow+eigvalhigh-1)>n) {
        Nume=n-eigvallow+1;
    } else Nume = eigvalhigh;
    SetDefault(Nume);

    iselec = eigval;   //new ntyp[limmax];
}

void Davidson::SetEigSelect(ntyp eigval[])
{
	if (!tab) delete[] iselec;
    iselec =  eigval;

    ilow = 0;
    ntyp eigvallow=eigval[0], eigvalhigh=eigval[0], Nume;
    for (ntyp i=1; eigval[i]>0; i++) {
    	eigvallow=MIN(eigvallow, eigval[i]);
        eigvalhigh=MAX(eigvalhigh, eigval[i]);
    }
    if ((eigvallow+eigvalhigh-1)>n) {
        Nume=n-eigvallow+1;
    } else Nume = eigvalhigh;
    SetDefault(Nume);

    return;
}

Davidson::~Davidson()
{
    if (!tab) delete[] iselec;
}

/* Default setup of davidson variables 
*/
void Davidson::SetDefault(ntyp Nume)
{
	Init=false;
    ierr=0;
    crite = 1e-15; critc = 1e-12; critr = 1e-8; ortho = 1e-9;
    niv = 0;
    mblock = 1;
    numemax = Nume;
    maxiter = MAX( numemax*40, 200 ); //maximum number of iterations by default Max(40*number_of_eig, or 200)
    nmv=0;
    neig = 0;

    limmax=numemax+20;
    iiwsz = 6*limmax + numemax;
    lim = limmax;

    irwsz = 2*n*limmax + limmax*limmax + (numemax+10)*limmax + numemax;

    hiend=false;

    norm=0;
    silence=false;
    loop=0;
    inf=0;
}

void Davidson::SetNumeMax(ntyp NumeMax)
{
    numemax=NumeMax;
    maxiter = MAX( numemax*40, 200 );
    ihigh = numemax;
    limmax=numemax+20;
    iiwsz = 6*limmax + numemax;
    lim = limmax;
    irwsz = 2*n*limmax + limmax*limmax + (numemax+10)*limmax + numemax;

    if (!tab) delete[] iselec;
    if (!tab) iselec = new ntyp[limmax];
}

ntyp Davidson::GetNumeMax() const
{
	return numemax;
}

ntyp Davidson::GetLimMax() const
{
	return limmax;
}

//Set the maximum limit for the size of 
void Davidson::SetLimMax(ntyp LimMax)
{
    limmax=LimMax;
    iiwsz = 6*limmax + numemax;
    lim = limmax;
    irwsz = 2*n*limmax + limmax*limmax + (numemax+10)*limmax + numemax;

    if (!tab) delete[] iselec;

    if (!tab) iselec = new ntyp[limmax];

}

void Davidson::CheckError(ntyp iselec[])
{
	ntyp err=0;
    if (lim>n)
        ++err;
    if (lim<=0) err += 2;

    if ((ilow<=0)||(ilow>n))
    {
    	//..Look for user choice of eigenpairs in ISELEC
    	if (iselec[0]<=0)
        	//..Nothing is given in ISELEC
        	err += 4;
        else {
        	/* ..Find number of eigenpairs wanted, and their            
               ..min/max indices    */
        	neig=1;
            ilow=iselec[0];
            ihigh=iselec[0];
            for (ntyp i=1; iselec[i]>0; i++) {
            	ilow=MIN(ilow, iselec[i]);
                ihigh=MAX(ihigh, iselec[i]);
                ++neig;
            }
            //..Check if a very large index is asked for
            if (ihigh>n) err += 8;
        }
    }
    else {
    	/* ..Look for a range between ILOW and IHIGH
           ..Invalid range. IHIGH>N */
    	if (ihigh>n) err += 8;
        neig=ihigh-ilow+1;
        //..Invalid range. IHIGH<ILOW
        if (neig<=0) err += 16;
        if (neig>lim)
        	//..Not enough Basis space. Increase LIM or decrease NEIG
        	err += 32;
        else
            //..Fill in the ISELEC with the required indices
        	for (ntyp i=1; i<=neig; i++) iselec[i-1]=ilow+i-1;
    }
    nume=ihigh;
    //..Identify if few of the highest eigenpairs are wanted.
    if ((ilow+ihigh-1)>n) {
    	hiend=true;
        nume=n-ilow+1;
    	/* ..Change the problem to a minimum eipenpairs one
           ..by picking the corresponding eigenpairs on the
           ..opposite side of the spectrum.    */
        for (ntyp i=1; i<=neig; i++) iselec[i-1]=n-iselec[i-1]+1;
    }
    if (neig>nume) err += 64;  //..duplications in ISELEC
    if ((nume>lim)||((nume==lim)&&(nume!=n))) err += 128; //..Not enough Basis space. Increase LIM or decrease NUME
    if ( (mblock<1)||(mblock>neig) ) err += 256; //..Size of Block out of bounds
    if ((irwsz<(lim*(2*n+lim+(nume+10))+nume)) ||
        (iiwsz<(6*lim+nume))) err += 512; //..Check for enough workspace for Dvdson

	if (niv>lim)
    	//..Check number of initial estimates NIV is lower than LIM.
    	cerr<<"\n\nWARNING: Too many initial estimates. \n "
        	  "The routine will pick the appropriate number \n";
    else if (niv<nume && niv>0)
    	//..check if enough initial estimates.
        //..(NIV<1 => program chooses)
    	cerr<<"\n\nWARNING: Not enough initial estimates \n"
        	  "The routine will pick the appropriate number \n";
}


void Davidson::Calculate()
{
    ExecuteDavidson();
}

//this is the main program 
void Davidson::ExecuteDavidson()
{
    if (!tab)
    	for (ntyp i=0; i<limmax; i++) iselec[i]=0;

    //set of initial vectors
    if (Init) InitVectors();
    
    CheckError(iselec);
    

    //Assigning space for the real work arrays (addressing)
//	ntyp ibasis    =0, //first we have lim vectors of size n
//         ieigval   =ibasis  +n*lim,			//then we have lim eigenvalues
//       	 iab       =ieigval +lim,			// then another set of lim vectors of size n (This is D=A*B, product set)
//     	 is        =iab     +n*lim,			//small matrix in krylov space, symmetric 
//    	 itemps    =is      +lim*(lim+1)/2, //scratch space, for diagonalization?
//    	 isvec     =itemps  +lim*(lim+1)/2, //array holding eigenvectors of S (only nume of them?) 
//    	 iscra1    =isvec   +lim*nume,       //scretch array used by DSPEVX
//    	 ioldval   =iscra1  +8*lim,     //Array keeping the previous' step eigenvalue estimates.
//		 
//    	 iscra2    =0, //Integer Srcatch array used by DSPEVX. //diagonalization
//    	 iscra3    =iscra2  +5*lim,
//    	 iicv      =iscra3  +lim;

    if (hiend) dscal(n,-1.0,diag->Diagonal(),1);
    ntyp istart=niv;
    
    //allocate matrices
    twosvec.resize(nume);ab.resize(lim);
    scratch1.resize(lim*nume,0);old.resize(lim,0);
    smt.resize(lim*(lim+1)/2);scra2.resize(5*lim,0);scra3.resize(lim,0);
    icv.resize(nume,0);
    for (int i=0;i<lim;i++)
    {
        ab[i].resize(n);
    }
    for (int i=0;i<nume;i++)
    {
        twosvec[i].resize(lim);
    }
    //allocate b if no vectors have been read in it
    if (b.size()==0)
    {
        b.resize(lim);
        eigs.resize(lim);
        for (int i=0;i<lim;i++)
            b[i].resize(n,0);
    }
    else
    {
        int ns=b.size();
        while(b.size()<lim)
        {
            b.push_back(ab[lim-1]);
            eigs.push_back(1.e+30);
        }

    }
    //add vectors to ab and b
    setup_2d(n, lim,nume,hiend,b,ab,smt,&istart,diag);
    nmv=istart;
    
    //start routine.
    dvdrvr2(n,hiend,lim,mblock, nume,istart,neig,iselec,crite,critc,critr,
           ortho,maxiter,eigs,b,ab,smt,twosvec,scratch1, scra2,scra3, icv,old, &nmv,&ierr, &loop, diag, silence, inf);
    if (hiend) {
    	dscal(n,-1.0,diag->Diagonal(),1);
        dscal(nume,-1.0,&eigs[0],1);
    }
/*  -Copy the eigenvalues after the eigenvectors
*   -Next, copy the difference of eigenvalues between the last two steps  
*   -Next, copy the residuals for the first NUME estimates */
    b[nume]=eigs;
    b[nume+1]=old;
    b[nume+2]=scratch1;
    return;
}

void Davidson::EigenValue(eig_val &Eigen, 
ntyp nx ///<eigenvalue number counting from zero
) const
{   
//    Eigen.Eigenvalue=work[nume*Davidson::n+nx];
//    Eigen.Eigval_differences=work[nume*Davidson::n+nume+nx];
//    Eigen.Residual=work[nume*Davidson::n+2*nume+nx];
    Eigen.Eigenvalue=b[nume][nx];
    Eigen.Eigval_differences=b[nume+1][nx];
    Eigen.Residual=b[nume+2][nx];
}

void Davidson::PrintEigenValue(ntyp nx) const
{   cout<<b[nume][nx]<<" "<<b[nume+1][nx]<<" "<<b[nume+2][nx]<<endl;
}

void Davidson::PrintEigenValue(void) const
{   
ntyp nvec=GetNumeMax();
for (ntyp i_nvec=0;i_nvec<nvec;i_nvec++) {PrintEigenValue(i_nvec);}
}
//save individual vector
void Davidson::SaveVector(std::ofstream &outputf, ntyp nx) const
{   
    uint_t_ee_file ui_n=(uint_t_ee_file) (Davidson::n);
    outputf.write((char*)(&(b[nume][nx])),sizeof(ftyp)); //save eigenvlaue
    outputf.write((char*)(&(b[nume+1][nx])),sizeof(ftyp));//save eigval err
    outputf.write((char*)(&(b[nume+2][nx])),sizeof(ftyp));//save eigvec err
	outputf.write((char*)(&ui_n),sizeof(uint_t_ee_file)); //save dimension
	outputf.write((char*)(&(b[(nx)][0])),(Davidson::n)*sizeof(ftyp)); //save eigenvector
}

void Davidson::Save(std::ofstream &outputf, ntyp nvec) const
{   
    unsigned ui_nvec=(unsigned) (nvec);
	outputf.write((char*)(&(ui_nvec)),sizeof(unsigned)); //save number of vectors
	for (ntyp i_nvec=0;i_nvec<nvec;i_nvec++) {SaveVector(outputf,i_nvec);}
}

void Davidson::Save(std::ofstream &outputf) const
{   
	ntyp nvec=GetNumeMax();
	Save(outputf, nvec);
}

void Davidson::Save(const char *outfile) const
{   
  ofstream outputf;
  outputf.open(outfile, ios::binary);
 if (!outputf) cerr<<"can not open output "<<outfile<<endl;
 else cerr<<" writing output to "<< outfile<<endl;
 Save(outputf);
 outputf.close();
}
//save some number of vectors
void Davidson::Save(const char *outfile, ntyp n) const
{   
  ofstream outputf;
  outputf.open(outfile, ios::binary);
 if (!outputf) cerr<<"can not open output "<<outfile<<endl;
 else cerr<<" writing output to "<< outfile<<endl;
 Save(outputf,n);
 outputf.close();
}

ntyp Davidson::GetIerr() const
{
	return ierr;//return the error flag
}

bool Davidson::GetHiend() const
{
	return hiend;// hiend is true if highest eigenpairs are sought
}

ntyp Davidson::GetMatrixVectProd() const
{
	return nmv;//return the number of matrix vector products.
}

//this code reads initial vectors
void Davidson::InitVectors()
{
    
    b.clear(); //vector of basis
    ifstream inputputf;
    inputputf.open(initfromfiles, ios::binary);
    uint_t_ee_file noinitvectors=0; //number of initial vectors
    inputputf.read((char*)(&(noinitvectors)),sizeof(unsigned)); //save number of vectors
    b.resize(noinitvectors); //vector of basis size by number
    eigs.resize(noinitvectors);
    cerr<<"Number of input vectors "<<noinitvectors<<endl;
//    read all vectors one by one
    ftyp EE,dE,dV; //temporary energy, error, variance
    uint_t_ee_file ndim; //read dimension
//    ftyp *p_work; //pointer to read array;
    for (int i=0; i<noinitvectors;i++) {
    b[i].resize(n); //create empty vector in basis
//        p_work=&(temp[0]);
        inputputf.read((char*)(&(EE)),sizeof(ftyp)); //read energy
        eigs[i]=EE;
        inputputf.read((char*)(&(dE)),sizeof(ftyp)); //read energy error
        inputputf.read((char*)(&(dV)),sizeof(ftyp)); //read energy variance
        cerr<<"State "<<i<<" Initial Energy   "<<EE<<" Errors:"<<dE<<" "<<dV<<endl;
        inputputf.read((char*)(&(ndim)),sizeof(uint_t_ee_file)); //read dimension
        if (ndim!=n) { cerr<<"Dimension mismatch matrix n="<<n<<" file n="<<ndim<<"\n"; exit (1);}
//        cerr<<"Read dimension "<<ndim<<endl;
        inputputf.read((char*)(&(b[i][0])),ndim*sizeof(ftyp)); //read the actual vector
//        p_work+=ndim; //shift the pointer position for next reading
//        temp.assign(n,0);
    }
    
    inputputf.close();
    niv=noinitvectors; //put the number of initial vectors as class variable
}

//we are now set to read initial vectors upon execution of davidson
void Davidson::ReadInitVectors(const char initfromfiles[])
{
	strcpy(Davidson::initfromfiles,initfromfiles);
    Init=true;
}


