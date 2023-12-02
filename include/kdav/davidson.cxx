//there were some trics involved, functions Multiply must be const (not modify internal class) for the
//virtual class to override 
#include <iostream>
#include "defaulttd.h"
#include "matrixbase.h"
#include "MandD.cxx" //diagonal to matrix
#include <cstdlib>
using std::atoi;
//#include "DavidsonBase.cxx"

#include "cdavidson.h"
#include "cdavidson.cpp"
#include "fsdavidson.cxx"




template <typename MyMatrixBase, typename evector, typename ematrix>
int DavidsonDiagonalizationBasic (
evector &EE, //<[out] eigenvalues, must have EE[i] access
ematrix &EV,//<[out] eigenvectors EV[i][0] must start contineous ftyp* memory block
const MyMatrixBase &FHH, //< matrix to diagonalize
int nstates=10, //number of states
int maxiter=0, //maximum number of iterations 0-set yourself
int savestates=0 //number of states to save
)
{

MyMatrixDiagonal<MyMatrixBase > *diag;
diag=new MyMatrixDiagonal<MyMatrixBase >(FHH);

char initvectors[255]="initvectors.par";
ntyp Init=0;

int select=0; //not set, eigen-pairs are computed for a set from 0 to nstates
//for other select we comput only particular states
ftyp crite=1e-15, critc=1e-12, critr=1e-8, ortho=1e-9;
ntyp numemax=10, limmax=0, mblock=10, norm=1;
//maxiter is maximum number of iterations, set below
numemax=mblock=nstates; //set block and number of states
Davidson *dav=0;
 eig_val myeig; //my eigenvalue
//ftyp EV[NN];
            if (select)
            	dav = new Davidson(diag, select);
            else
        		dav = new Davidson(diag,numemax,1);
            //this allows to overwrite the variables set by default in cdavidson.cpp::void Davidson::SetDefault(ntyp Nume) 
            dav->crite=crite; dav->critc=critc; dav->critr=critr; dav->ortho=ortho;
			dav->silence=0;
            if (maxiter==-1) dav->inf=true;
            else if (!maxiter) maxiter=MAX(dav->GetNumeMax()*40, 200);
			//GetNumeMax() tells the highest eigen-pair to be computed
            dav->maxiter=maxiter;
            if (!limmax) limmax=dav->GetNumeMax()+20;
            dav->SetLimMax(limmax); //maximum basis size
            dav->mblock=mblock;
//			cerr<<endl;
//			cerr<<"Davidson: Block: "<<numemax<<" Iterations: "<<maxiter<<" Eigenvalues: ";
//			if (savestates!=0) cerr<<savestates<<endl; else cerr<<numemax<<endl;
            if (Init) dav->ReadInitVectors(initvectors);
            dav->Calculate();
//			if (savestates!=0) dav->Save(outf, savestates); //save desired number of states
//			else dav->Save(outf);
			//dav->EigenValue(myeig,0);
//			dav->PrintEigenValue();
            for (int s=0;s<savestates;s++) {
            dav->EigenValue(myeig,s); //get data in eigenvalue
            EE[s]=myeig.Eigenvalue; //save eigenvalue
//                cout<<EE[0]<<endl;
//			dav->EigenFunction(s,&(EV[s][0]),1); //save normalized eigenvector
                //eigenvectors come out normalised to 1.
            for (int i=0;i<(dav->n);i++)
                EV[s][i]=dav->b[s][i];
            }
  //  delete dav;
  //  delete diag; should not have delete function
			//dav->EigenFunction(2,EV,1);
					
		//cout<<"my lowest eigenvalue is "<<myeig.Eigenvalue<<endl;

return 0;		

		}


/*! \brief This also reads initial vectors must be the same number as nstates
 */
template <typename MyMatrixBase, typename evector, typename ematrix>
int DavidsonDiagonalizationBasic (
                                  evector &EE, //<[out] eigenvalues, must have EE[i] access
                                  ematrix &EV,//<[out] eigenvectors EV[i][0] must start contineous ftyp* memory block
                                  ematrix &INEV,//<[in] eigenvectors IN[i][0] must start contineous ftyp* memory block
                                  const MyMatrixBase &FHH, //< matrix to diagonalize
                                  int nstates=10, //number of states
                                  int maxiter=0, //maximum number of iterations 0-set yourself
                                  int savestates=0 //number of states to save
                                  )
{
    
    MyMatrixDiagonal<MyMatrixBase > *diag;
    diag=new MyMatrixDiagonal<MyMatrixBase >(FHH);
    
    //char initvectors[255]="initvectors.par";
    ntyp Init=1; //set read initial vectors
        
    int select=0; //not set, eigen-pairs are computed for a set from 0 to nstates
    //for other select we comput only particular states
    ftyp crite=1e-15, critc=1e-12, critr=1e-8, ortho=1e-9;
    ntyp numemax=10, limmax=0, mblock=10, norm=1;
    //maxiter is maximum number of iterations, set below
    numemax=mblock=nstates; //set block and number of states
    Davidson *dav=0;
    eig_val myeig; //my eigenvalue
    //ftyp EV[NN];
    if (select)
        dav = new Davidson(diag, select);
    else
        dav = new Davidson(diag,numemax,1);
    //this allows to overwrite the variables set by default in cdavidson.cpp::void Davidson::SetDefault(ntyp Nume)
    dav->crite=crite; dav->critc=critc; dav->critr=critr; dav->ortho=ortho;
    dav->silence=0;
    if (maxiter==-1) dav->inf=true;
    else if (!maxiter) maxiter=MAX(dav->GetNumeMax()*40, 200);
    //GetNumeMax() tells the highest eigen-pair to be computed
    dav->maxiter=maxiter;
    if (!limmax) limmax=dav->GetNumeMax()+20;
    dav->SetLimMax(limmax); //maximum basis size
    dav->mblock=mblock;
    if (Init) { //dav->ReadInitVectors(initvectors);
        dav->Init2=true;
        dav->b.resize(nstates);
        
        dav->niv=nstates; //put the number of initial vectors as class variable
        for (int i=0;i<nstates;i++)
        {
            dav->b[i].resize(dav->n);
            for (ntyp q=0;q<(dav->n);q++)
            {
                dav->b[i][q]=INEV[i][q];
            }
        }
    }
    dav->Calculate();
    //			if (savestates!=0) dav->Save(outf, savestates); //save desired number of states
    //			else dav->Save(outf);
    //dav->EigenValue(myeig,0);
    //			dav->PrintEigenValue();
    for (int s=0;s<savestates;s++) {
        dav->EigenValue(myeig,s); //get data in eigenvalue
        EE[s]=myeig.Eigenvalue; //save eigenvalue
//        dav->EigenFunction(s,&(EV[s][0]),1); //save normalized eigenvector
        for (int i=0;i<(dav->n);i++)
            EV[s][i]=dav->b[s][i];
    }
    
    //dav->EigenFunction(2,EV,1);
	
    //cout<<"my lowest eigenvalue is "<<myeig.Eigenvalue<<endl;
    
    return 0;		
}


/*! \brief Templated matrix diagonalizaiton with output file
*/
template <typename MyMatrixBase>
int DavidsonDiagonalization (
const char * outf, //<output file name
const MyMatrixBase &FHH, //< matrix to diagonalize
int nstates=10, //number of states
int maxiter=0, //maximum number of iterations 0-set yourself
int savestates=0 //number of states to save
)
{

MyMatrixDiagonal<MyMatrixBase > *diag;
diag=new MyMatrixDiagonal<MyMatrixBase >(FHH);

char initvectors[255]="initvectors.par";
ntyp Init=0;




int select=0; //not set, eigen-pairs are computed for a set from 0 to nstates
//for other select we comput only particular states
ftyp crite=1e-15, critc=1e-12, critr=1e-8, ortho=1e-9;
ntyp numemax=10, limmax=0, mblock=10, norm=1;
//maxiter is maximum number of iterations, set below
numemax=mblock=nstates; //set block and number of states
Davidson *dav=0;
 eig_val myeig; //my eigenvalue
//ftyp EV[NN];
    cerr<<"Size: "<<diag->Size()<<endl;
            if (select)
            	dav = new Davidson(diag, select);
            else
        		dav = new Davidson(diag,numemax,1);
            //this allows to overwrite the variables set by default in cdavidson.cpp::void Davidson::SetDefault(ntyp Nume) 
            dav->crite=crite; dav->critc=critc; dav->critr=critr; dav->ortho=ortho;
			dav->silence=0;
            if (maxiter==-1) dav->inf=true;
            else if (!maxiter) maxiter=MAX(dav->GetNumeMax()*40, 200);
			//GetNumeMax() tells the highest eigen-pair to be computed
            dav->maxiter=maxiter;
            if (!limmax) limmax=dav->GetNumeMax()+20;
            dav->SetLimMax(limmax); //maximum basis size
            dav->mblock=mblock;
			cerr<<endl;
			cerr<<"Davidson: Block: "<<numemax<<" Iterations: "<<maxiter<<" Eigenvalues: ";
			if (savestates!=0) cerr<<savestates<<endl; else cerr<<numemax<<endl;
            if (Init) dav->ReadInitVectors(initvectors);
            dav->Calculate();

			if (savestates!=0) dav->Save(outf, savestates); //save desired number of states
			else dav->Save(outf);
			dav->EigenValue(myeig,0);
			dav->PrintEigenValue();
		cerr<<"my lowest eigenvalue is "<<myeig.Eigenvalue<<endl;

return 0;		
		}

//same as above but reads initial vectors
template <typename MyMatrixBase>
int DavidsonDiagonalization (
                             const char * outf, //<output file name
                             const MyMatrixBase &FHH, //< matrix to diagonalize
                             const char* initvectors, //filename that stores initial vectors in .EE
                             int nstates, //number of states
                             int maxiter, //maximum number of iterations 0-set yourself
                             int savestates //number of states to save
                             )
{
    
    MyMatrixDiagonal<MyMatrixBase > *diag;
    diag=new MyMatrixDiagonal<MyMatrixBase >(FHH);
    
    //char initvectors[255]="initvectors.par";
    ntyp Init=1;
    
    int select=0; //not set, eigen-pairs are computed for a set from 0 to nstates
    //for other select we comput only particular states
    ftyp crite=1e-15, critc=1e-12, critr=1e-8, ortho=1e-9;
    ntyp numemax=10, limmax=0, mblock=10, norm=1;
    //maxiter is maximum number of iterations, set below
    numemax=mblock=nstates; //set block and number of states
    Davidson *dav=0;
    eig_val myeig; //my eigenvalue
    //ftyp EV[NN];
        if (select)
            dav = new Davidson(diag, select);
        else
            dav = new Davidson(diag,numemax,1);
        //this allows to overwrite the variables set by default in cdavidson.cpp::void Davidson::SetDefault(ntyp Nume) 
        dav->crite=crite; dav->critc=critc; dav->critr=critr; dav->ortho=ortho;
        dav->silence=0;
        if (maxiter==-1) dav->inf=true;
        else if (!maxiter) maxiter=MAX(dav->GetNumeMax()*40, 200);
        //GetNumeMax() tells the highest eigen-pair to be computed
        dav->maxiter=maxiter;
        if (!limmax) limmax=dav->GetNumeMax()+20;
        dav->SetLimMax(limmax); //maximum basis size
        dav->mblock=mblock;
        cerr<<endl;
        cerr<<"Davidson: Block: "<<numemax<<" Iterations: "<<maxiter<<" Eigenvalues: ";
        if (savestates!=0) cerr<<savestates<<endl; else cerr<<numemax<<endl;
        if (Init) dav->ReadInitVectors(initvectors);
        dav->Calculate();
        if (savestates!=0) dav->Save(outf, savestates); //save desired number of states
        else dav->Save(outf);
        dav->EigenValue(myeig,0);
        dav->PrintEigenValue();
        //dav->EigenFunction(2,EV,1);

    
    cerr<<"my lowest eigenvalue is "<<myeig.Eigenvalue<<endl;

    return 0;		
}
