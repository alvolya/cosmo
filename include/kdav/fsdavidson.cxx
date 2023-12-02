//#define STATUS yes
//#include <status.h>
#include "defaulttd.h"
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>
using namespace std;
#include <tav/EigenSystem.cxx>
//#include <tav/sort.cxx>
#include <algorithm>
#include "matrixbase.h"

#include <fstream>
#ifndef  __FSDAVIDSON_CXX
#define __FSDAVIDSON_CXX

#define MAX(x,y) (((x)>(y))? (x) : (y))
#define MIN(x,y) (((x)<(y))? (x) : (y))
//---------------------------------------------------------------------------
void DiagonalizeS(vector<ftyp> &s,vector<vector<ftyp> > &twosvec,vector<ftyp> &eigval,ntyp basis_dim,ntyp nume)
{

    vector<ntyp> index(basis_dim,0);
    twosvec.resize(basis_dim);
    for (ntyp i=0;i<basis_dim;i++)
    {
        index[i]=i;

        twosvec[i].resize(basis_dim,0);
    }

    for (ntyp i=0;i<basis_dim;i++)
    {
        for (ntyp j=0;j<=i;j++)
        {
            twosvec[i][j]=s[i*(i+1)/2+j];
            twosvec[j][i]=s[i*(i+1)/2+j];
        }
    }

    tav::EigenSystemT(basis_dim, twosvec, eigval);
    //tav::QuickSort(basis_dim,eigval,index);
    // Sort the index array based on the values in eigval
    //index is filled before
    //     for (int i = 0; i < basis_dim; i++) {
    //     index[i] = i;
    // }
    //use [&eigval] if vector and just [eigval] if array pointer
    std::sort(index.begin(), index.end(), [&eigval](int i1, int i2) {
       return eigval[i1] < eigval[i2];
    });
    //sort(eigval, eigval + basis_dim); //things are too small so we can sort eigval again keeping index
    std::sort(eigval.begin(), eigval.end());
    
    // Create a temporary array to hold the sorted values
    // vector<double> sortedEigval(basis_dim);
    // for (int i = 0; i < basis_dim; i++) {
    //     sortedEigval[i] = eigval[index[i]];
    // }

    // // Copy the sorted values back to eigval
    // for (int i = 0; i < basis_dim; i++) {
    //     eigval[i] = sortedEigval[i];
    // }

// //sort using algorithm function---------------------------------
//  // Pair each element of eigval[] with its index
//     vector<pair<double, ntyp>> paired(basis_dim);
//     for (auto i = 0; i < basis_dim; i++) {
//         paired[i] = {eigval[i], i};
//     }

//     // Sort the pair array in ascending order
//     sort(paired.begin(), paired.end());

//     // Update the original arrays with sorted values and indices
//     for (auto i = 0; i < basis_dim; i++) {
//         eigval[i] = paired[i].first;
//         index[i] = paired[i].second;
//     }
//     //this replaces the quicksort used above-----------------

    for (ntyp i=0;i<nume;i++)
    {
        twosvec[i].swap(twosvec[index[i]]);
        for (ntyp j=i;j<nume;j++)
        {
            if (index[j]==i)
                index[j]=index[i];
        }
    }
    ntyp si=twosvec.size();
    for (int i=nume;i<si;i++)
    {
        twosvec.pop_back();
    }
}
void Collapse_Matrix(  vector<vector<ftyp> > &c,vector<vector<ftyp> > &b,vector<double*> &pt,int kp)
{
    ntyp n=c.size(),N=b[0].size();
#pragma omp parallel for schedule(static)
    for (ntyp i=0;i<N;i++)
    {
        vector<double> temp(n,0);
        for (int j=0;j<n;j++)
        {
            temp[j]=0;
            for (int k=0;k<kp;k++)
            {
                temp[j]+=c[j][k]*(*(pt[k]+i));//b[k][i];
            }
//            *(pt[j]+i)=temp[j];

        }
        for (int j=0;j<n;j++)
            *(pt[j]+i)=temp[j];
        
    }
    return;
}
double dot_vectors(vector<double> &dx,  vector<double> &dy)
{
    ftyp dtemp=0.0;
    ntyp n=dy.size();
//    for (int i=0;i<n;i++)
//        dtemp+=dx[i]*dy[i];
//    return dtemp;
    ntyp m=n%5;
    if (m){
        for (ntyp i=1; i<=m; i++) dtemp += dx[i -1]*dy[i -1];
        if (n<5) return dtemp;
    }
    for (ntyp i=m + 1; i<=n; i+=5)
        dtemp += dx[i -1]*dy[i -1] + dx[i]*dy[i] + dx[i + 1]*dy[i + 1]
        + dx[i + 2]*dy[i + 2] + dx[i + 3]*dy[i + 3];
    return dtemp;
}
void scale_vector(vector<double> &a,double b,ntyp N)
{
#pragma omp parallel for schedule(static)
    for(ntyp i=0;i<N;i++)
        a[i]*=b;
}
ntyp Find_Max_Element_Index( ntyp n///<size of array
                                        ,ftyp dx[]///<the array
                                        )
{
    ftyp dmax=fabs(dx[0]);
    ntyp ix=0,fidamax=1;
    for (ntyp i=2; i<=n; i++)
    {
      	if (fabs(dx[i -1])>dmax)
        {
        	fidamax=i;
            dmax=fabs(dx[i -1]);
        }
    }
    return fidamax;
}
void Add_Vectors_To_AB_S(ntyp n, ntyp lim, bool hiend, ntyp basis_dim,
                                     ntyp nncv, //<new dimension for the dot multiplication
                                     vector<vector<ftyp> > &basis, vector<vector<ftyp> > &ab, vector<ftyp>  &s,
                                     const MatB::MatrixBase *const diag)
{
    // The user specified matrix-vector routine is called with the new
    // basis vector B(*,basis_dim+1) and the result is assigned to AB(idstart)
    //    ntyp idstart=kpass*n;
    ntyp idstart=basis_dim;
    // diag->dssbmv(nncv,&basis[idstart -1],&ab[idstart -1],n);
    //	diag->TimesMatrix_2d(ab,basis,nncv); //new multiplication no need for n-dimension of the matrix
    
    double **y=new double* [nncv];
    double **x=new double* [nncv];
    for (ntyp i=idstart;i<(idstart+nncv);i++)
    {
        y[i-idstart]=&(ab[i][0]);
        x[i-idstart]=&(basis[i][0]);
    }
    
    diag->TimesMatrix(y,x,nncv);
    
    //    // If highest pairs are sought, use the negative of the matrix
    if (hiend)
    {
        for (ntyp i=0;i<nncv;i++)
            scale_vector(ab[basis_dim+i],-1.0,ab[0].size());
    }
    //
    //    // The new S is calculated by adding the new last columns
    //    // S(new)=B^T D(new).
    ntyp isstart=basis_dim*(basis_dim+1)/2;
    ntyp iv=0,ibv=0;
    for (iv=1; iv<=nncv; iv++)
    {
    	ntyp ibstart=0;
#pragma omp parallel for schedule(static)
        for (ibv=0; ibv<basis_dim+iv; ibv++)
        {
            //           	s[isstart + ibv]=dot_vectors(basis[ibstart],ab[idstart]);
            s[isstart + ibv]=dot_vectors(basis[ibv],ab[idstart]);
            //           	ibstart++;
        }
        isstart += (basis_dim+iv);
        idstart++;
    }
    
    delete [] y;
    delete [] x;
	return;
}
bool Test_Convergence( ntyp basis_dim///<Current dimension of basis
                                  ,ntyp nume
                                  ,ntyp neig
                                  ,vector<ntyp> &iselec///<indices of sought eigenvectors
                                  ,vector<vector<ftyp> > &svec///<Eigenvectors of S=B^T A B
                                  ,vector<ftyp> &eigval///<Eigenvalues
                                  ,vector<ntyp> &icv///<icv[i]=0 if i-th eigenpair has converged, 0 if not
                                  ,ftyp crite///<Convergenvce Threshold for eigenvalues
                                  ,ftyp critc///<Convergence Threshold for eigenvectors
                                  ,vector<double> &oldval///<Eigenvalues from previous iteration
                                  ,ntyp *nncv///<Number of non-converged eigenpairs
                                  ,vector<ntyp> &incv///<indices of sought eigenpairs in order of decreasing maximum coefficient
                                  )

{
    bool done=false;
    vector<double> rowlast(neig,0);
    vector<int> ind(neig,0);
    // Test all wanted eigenvalues for convergence under CRITE
    ntyp nnce=0;
    for (ntyp i=1; i<=neig; i++) {
    	ntyp ival=iselec[i -1];
        if (fabs(oldval[ival -1]-eigval[ival -1])>=crite) ++nnce;
    }
    if (!nnce) return true;
    
    // Find the maximum element of the last NNCV coefficients of unconverged
    // eigenvectors. For those unconverged coefficients, put their indices
    // to IND and find their number NNCV
    ftyp tmax=0, temp=0;
    ntyp indx=0, itemp=0, icur=0, icnt=0;
    for (ntyp i=1; i<=neig; i++)
    	if (icv[iselec[i -1] -1]==0)
        {
        	//..Find coefficient and test for convergence
        	icur=basis_dim*iselec[i -1];
            tmax=fabs(svec[iselec[i-1]-1][basis_dim-1]);
            for (ntyp l=1; l<=*nncv-1; l++)
            {
                tmax=MAX(tmax, fabs(svec[iselec[i-1]-1][basis_dim-l -1]));
            }
            if (tmax<critc)//..this  coefficient converged
            {
                icv[iselec[i -1] -1]=1;
            }
            else
            {
            	//..Not converged. Add it to the list.
            	++icnt;
                ind[icnt -1]=iselec[i -1];
                rowlast[icnt -1]=tmax;
            }
        }
    *nncv=icnt;
    if (!(*nncv)) done=true;
    
    // Sort the ROWLAST elements interchanging their indices as well
    for (ntyp i=1; i<=*nncv; i++) {
    	indx=Find_Max_Element_Index(*nncv-i+1,&rowlast[i -1]);
        incv[i -1]=ind[indx+i-1 -1];
        temp=rowlast[indx+i-1  -1];
        rowlast[indx+i-1 -1]=rowlast[i -1];
        rowlast[i -1]=temp;
        itemp=ind[indx+i-1 -1];
        ind[indx+i-1 -1]=ind[i -1];
        ind[i -1]=itemp;
    }
    
	return done;
}
void daxpy_vector(ftyp da,vector<ftyp> &dx,vector<ftyp> &dy)
{
    // y=y+a*x (for vectors y,x and scalar a)
    ntyp n=dy.size();
    for (int i=0;i<n;i++)
        dy[i]+=da*dx[i];
//    ntyp m=n%4;
//    if (m){
//        for (ntyp i=1; i<=m; i++) dy[i-1]+=da*dx[i-1];
//        if (n<4) return;
//    }
//#pragma omp parallel for schedule(static)
//    for (ntyp i=m + 1; i<=n; i+=4){
//        dy[i-1] += da*dx[i-1];
//        dy[i] += da*dx[i];
//        dy[i + 1] += da*dx[i + 1];
//        dy[i + 2] += da*dx[i + 2];
//    }
}
//void The_Full_Daxpy(vector<double> &scra1,vector<vector<double> > &basis,int icur,int limit)
//{
//    ntyp i=0;
//    ntyp N=basis[icur].size();
//#pragma omp parallel for schedule(static)
//    for (i=0;i<N;i++)
//    {
//        for (int j=0;j<limit;j++)
//            basis[icur][i]+=-scra1[j]*basis[j][i];
//    }
//    
//}
void Orthogonalize_Basis(ntyp lim///<Limit to size of basis
                                     , ftyp ortho///<Limit below which orthogonality is broken
                                     , ntyp basis_dim///<Current dimension of basis
                                     , ntyp *nncv///<Number of new basis vectors
                                     , vector<ftyp> &scra1///<Scratch vector
                                     , vector<vector<ftyp> > &basis///<The basis Matrix
                                     , bool *restart///<If true the basis will be collapsed on return
                                     )
{
	bool powrot=false;
    ntyp  icur=0, iv=0,k=0;
    ftyp dprev=0, dcur=0;
    
    // ORTHOGONALIZATION
	*restart=false;
    icur=basis_dim;//set icur to first new vector.
    //.. do iv=1,nncv
    iv = 1;//loop over the number of newly added vector.
    do {
    	dprev=1.0e+7;//the maximum overalp of the previous vector.
        do {
        	powrot=false;
            dcur=0.0;
            
            //Take inner product of the new basis vector with all previous ones and store
            //the result in scratch space.
//            scra1.assign(scra1.size(),0);
#pragma omp parallel for schedule(static)
            for (k=1; k<=basis_dim+iv-1; k++)
            {
            	scra1[k-1]=dot_vectors(basis[k-1],basis[icur]);
            }
            dcur=fabs(scra1[0]);
            for (ntyp i=1;i<=basis_dim+iv-1;i++)
            {
                dcur=MAX(dcur,fabs(scra1[i-1]));//Store the maximum overlap between the new vector and all the old ones.
            }
            //Store the new orthogonal vector (icur)
//            The_Full_Daxpy(scra1,basis,icur,basis_dim+iv);
            for (ntyp i=1; i<=basis_dim+iv-1; i++)
            {
             	daxpy_vector(-scra1[i-1],basis[i-1],basis[icur]);
            }
            if (dcur>=ortho)//if the loss of precision is too great, restart.
            {
            	if (dcur>dprev)
                {
                	*restart=true;
                    //..Adjust the number of added vectors.
                    *nncv=iv-1;
                    return;
                }
                else
                {
                	dprev=dcur;
                    powrot=true;
                }
            }
        } while (powrot);
        
        // NORMALIZATION
        scra1[0]=dot_vectors(basis[icur],basis[icur]);//get the norm of the vector.
        scra1[0]=sqrt(scra1[0]);
        if (scra1[0]<1e-14) //if the norm is very small (vector is close to 0)
        {
            basis[icur]=basis[*nncv-1];//make it equal to last of the added vectors
            --*nncv;//decrease the number of added vectors by one.
        }
        else
        {
        	scale_vector(basis[icur],1./double(scra1[0]),basis[0].size());//rescale the vector by dividing with the norm
            ++icur;//move on to next vector
          	++iv;//include the vector just orthogonalised for the next one.
        }
    } while (iv<=*nncv);
    
    return;
}
void Matrix_Times_Vector_Plus_Vector(  vector<vector<ftyp> > &A,vector<ftyp> &x,vector<ftyp> &y,  ftyp alpha,  ftyp beta,ntyp dimension)
/*------------------------------------------------------------------------
 * called by: newvec_2d,
 *  same as dgemv but with vectors, for x,y and vector<vector> for A.
 *  if q=n normal matrix multiply with vector +vector
 *  a * A * x + b * y.
 * The result is stored in y.
 * if q=t same thing but transpose matrix A.
 *-----------------------------------------------------------------------*/
{
    if (beta!=0)
    {
#pragma omp parallel for schedule(static)
        for (int i=0;i<y.size();i++)
            y[i]*=beta;
    }
    else
    {
#pragma omp parallel for schedule(static)
        for (int i=0;i<y.size();i++)
            y[i]=0.0;
    }
#pragma omp parallel for schedule(static)
    for (ntyp i=0;i<A[0].size();i++)
    {
        ftyp *xx=NULL;
        xx=&x[0];
        
        for (ntyp j=0;j<dimension;j++)
        {
            y[i]+=A[j][i]*(*xx);
            xx++;
        }
    }
}
void Add_Vectors_To_Basis( ntyp nume///<index of highest eigenpair sought
                                      ,ntyp lim///<Limit to the basis size
                                      ,ntyp mblock///<Number of Targeted vectors
                                      ,ntyp basis_dim///<Current Dimension of basis
                                      ,ftyp critr///<Vector convergence threshold
                                      ,ftyp ortho///<orthogonality threshold
                                      ,ntyp *nncv///<number of unconverged eigenpairs
                                      ,vector<ntyp> &incv//indices of vectors to be targeted
                                      ,ftyp diag[]///<diagonal elements of Original Matrix for preconditioning
                                      ,vector<vector<ftyp> > &svec///<Eigenvectors of S=B^T A B
                                      ,vector<ftyp> &eigval///<eigenvalues of S=B^T A B
                                      ,vector<vector<ftyp> > &ab///<Product of Matrix A with the basis
                                      ,vector<vector<ftyp> > &basis///<The Basis Matrix
                                      ,vector<ntyp> &icv///<Index to wheter a vector has converged or not
                                      ,bool *restart
                                      ,bool *done)
{
	ftyp ss=0;
    ntyp indx=0,n=basis[0].size();
	*done   = false;
    ntyp nadded  = 0,
    icvc    = 0,
    limadd= MIN( lim, mblock+basis_dim ),/*In case adding mblock takes us over the maximum number of basis dimension allowed*/
    icur=basis_dim;//This is the point where the new vectors should start being placed
    // Compute RESIDUALS for the MBLOCK of the NNCV not converged vectors.
    //    gettimeofday(&start,NULL);
    for (ntyp i=0; i<*nncv; i++) {
    	indx=incv[i];
        //        indx=i;
        //         ..Compute  Newv=BASIS*Svec_indx , then                       ACPZ0893
        //         ..Compute  Newv=AB*Svec_indx - eigval*Newv and then          ACPZ0894
        //         ..compute the norm of the residual of Newv
        
        //Perform B^T tempor and store the result in B[icur]
        Matrix_Times_Vector_Plus_Vector(basis,svec[indx-1],basis[icur],1.0,0.0,basis_dim);
        Matrix_Times_Vector_Plus_Vector(ab,svec[indx-1],basis[icur],1.0,-eigval[indx-1],basis_dim);
        //Perform (AB)^T tempor - lambda[indx-1]*B[icur] and store it in B[icur]
        //Get the norm of the newly created vector
        ss = sqrt(dot_vectors(basis[icur],basis[icur]));
        
        //        //..Check for convergence of this residual
//        if ((indx)>=10)
//        {
//            cout<<indx-1<<" "<<ss<<" "<<eigval[indx-1]<<endl;
//        }
        if (ss<critr)
        {
//            cout<<indx-1<<" "<<ss<<endl;
            //        	//..Converged,do not add. Go for next non converged one
        	++icvc;
          	icv[indx-1] = 1;//Set icv to 1 => [indx-1] vector has converged.
            if (icvc<*nncv) continue;
            //            //..All have converged.
            *done=true;
          	return;
        }
        else
        {
        	//..Not converged. Add it in the basis
        	++nadded;
            incv[nadded -1]=indx;
            if ((nadded+basis_dim)==limadd) break;
            //..More to be added in the block
            icur++;//Move on to the next vector to add.
        }
    }
    *nncv=nadded;
    // Diagonal preconditioning: newvect(i)=newvect(i)/(l-Aii)               ACPZ0924
    // If (l-Aii) is very small (or zero) divide by 10.D-6
    icur=basis_dim;//move icur back to original place.
    ntyp i,irow;
//    ftyp dg;
    for (i=0; i<*nncv; i++)//For all new vectors
    {
#pragma omp parallel for schedule(static)
    	for (irow=0; irow<n; irow++) //for all elements of the vector
        {
            ftyp dg=eigval[incv[i]-1]-diag[irow];//substract the diagonal element of the Original Matrix (A) from the eigenvalue associated with the specific (i-th) new vector.
            //precondition.
            if (fabs(dg)>(1.0e-13))
                basis[basis_dim+i][irow] = basis[basis_dim+i][irow] / dg;
            else
                basis[basis_dim+i][irow]=basis[basis_dim+i][irow] /1.0e-13;
        }
    }
    
    // Orthonormalize the new vectors to the already existing basis.
    Orthogonalize_Basis(lim,ortho,basis_dim,nncv,ab[basis_dim], basis,restart);
//    cout<<"Restart: "<<*restart<<endl;
    return;
}
void Collapse_S(vector<vector<double> > &twosvec,vector<double> &stm,vector<double> &eigs,int nume, int lim)
{
    for (int i=0;i<lim*(lim+1)/2;i++)
        stm[i]=0.;
    ntyp ind=0;
    for (int i=1;i<=nume;i++)
    {
        stm[ind+i-1]=eigs[i-1];
        ind+=i;
    }
    for (int i=0;i<nume;i++)
    {
        
        for (int j=0;j<nume;j++)
        {
            if (i==j)
                twosvec[i][j]=1.;
            else
                twosvec[i][j]=0.;
            
        }
    }
}
void setup_2d(ntyp n, ntyp lim, ntyp nume, bool hiend, vector<vector<ftyp> > &basis, vector<vector<ftyp> > &ab, vector<ftyp> &s, ntyp *niv,const MatB::MatrixBase *const diag)
/*       Subroutine for setting up (i) the initial BASIS if not provided,
 *       (ii) the product of the matrix A with the Basis into matrix AB,
 *       and (iii) the small matrix S=B^TAB. If no initial estimates are
 *       available, the BASIS =(e_i1,e_i2,...,e_iNUME), where i1,i2,...,
 *       iNUME are the indices of the NUME lowest diagonal elements, and
 *       e_i the i-th unit vector. (ii) and (iii) are handled by ADDABS.
 *
 *
 *   on entry
 *   --------
 *   OP          The block matrix-vector operation, passed to ADDABS
 *   N           the order of the matrix A
 *   LIM         The limit on the size of the expanding Basis
 *   NUME        Largest index of the wanted eigenvalues.
 *   HIEND       Logical. True only if the highest eigenpairs are needed.
 *   diag        poiter for object of MatrixBase class
 *   MINELEM     Array keeping the indices of the NUME lowest diagonals.
 *
 *   on exit
 *   -------
 *   BASIS       The starting basis.
 *   AB, S       The starting D=AB, and small matrix S=B^TAB
 *   NIV         The starting dimension of BASIS.
 */
{
//    cout<<(*niv)<<" "<<lim<<" "<<nume<<endl;
	if ((*niv>lim) || (*niv<nume)) {
        vector<ntyp> minelem(nume,0);
        /*         ..Initial estimates are not available. Give as estimates unit
         *          ..vectors corresponding to the NUME minimum diagonal elements
         *          ..First find the indices of these NUME elements (in MINELEM).
         *          ..Array AB is used temporarily as a scratch array.	*/
//    	dinit(n,-1.0,ab,1);
        for (int i=0;i<n;i++)
            ab[0][i]=-1.;
        for (ntyp i=1; i<=nume; i++) {
        	ntyp imin=n+1;
            for (ntyp j=1; j<=n; j++)
            	if (ab[0][j -1]<0) {
                	imin=j;
                    break;
                }
            for (ntyp j=imin+1; j<=n; j++)
            	if ((ab[0][j-1]<0)&&((*diag)[j -1]<(*diag)[imin-1])) imin=j;
            minelem[i-1]=imin;
          	ab[0][imin-1]=1.0;
        }
        // ..Build the Basis. B_i=e_(MINELEM(i)) changed to C indexing
//        dinit(n*lim,0.0,basis,1);

        for (ntyp j=(*niv); j<nume; j++) {
        	ntyp i=ntyp((j)*n+minelem[j]);
            basis[j][minelem[j]-1]=1.0;
//            cout<<minelem[j]-1<<endl;
        }
        if ((*niv)!=0)
        {
            bool restart=false;
            ntyp a=nume-(*niv);
            Orthogonalize_Basis(basis.size(),1.e-9,(*niv),&a,ab[nume+1], basis,&restart);
        }
        *niv=nume;
    }
    
    /* Find the matrix AB by matrix-vector multiplies, as well as the        ACPZ0403
     *  small matrix S = B^TAB.	*/
    ntyp kpass=0;
    Add_Vectors_To_AB_S(n,lim,hiend,kpass,*niv,basis,ab,s,diag);
    
    return;
}



void dvdrvr2(ntyp n, bool hiend, ntyp lim, ntyp mblock, ntyp nume, ntyp niv, ntyp neig,
                        ntyp iselec[], ftyp crite, ftyp critc, ftyp critr, ftyp ortho, ntyp maxiter,
                        vector<ftyp> &eigs, vector<vector<ftyp> > &b, vector<vector<ftyp> > &ab, vector<ftyp > &stm,
                        vector<vector<ftyp> > &twosvec, vector<ftyp> &scra1, vector<ntyp> &iscra2, vector<ntyp> &incv2, vector<ntyp> &icv2, vector<ftyp> &old,
                        ntyp *nmv, ntyp *ierr, ntyp *loop,
                        const MatB::MatrixBase *const diag, bool silence, bool inf)
{
	bool restart=false, first=true, done=false, omin_to=false, userbreak=false;
    ntyp nfound=0, info=0, nloops=0;;
    ntyp ind=0;
    vector<ntyp> isel;
    vector<double*> ptb(b.size(),NULL),ptab(ab.size(),NULL);
    for (int i=0;i<b.size();i++)
    {
        ptb[i]=&b[i][0];
        ptab[i]=&ab[i][0];
    }
//    cout<<niv<<endl;
    ntyp kpass=niv,nncv=kpass;
    for (int i=0;i<neig;i++)
        isel.push_back(iselec[i]);
//    cout<<kpass<<" "<<nncv<<" "<<nume<<" "<<mblock<<endl;
#ifdef CSTATUS
    av::status S_D(maxiter,"D");
#ifdef FSTATUS
    S_D.SetFile();
#endif
#endif
    do {
#ifdef CSTATUS
        S_D(nloops);
#endif

//        cout<<kpass<<" "<<nncv<<" "<<nume<<endl;

        /*      (iterations for kpass=NUME,LIM)
         *
         * 		Diagonalize the matrix S. Find only the NUME smallest eigenpairs 	*/
//    	dcopy(nume,eigval,1,oldval,1);
//        cout<<"here0"<<endl;getchar();

        for (int i=0;i<nume;i++)
        old[i]=eigs[i];
////        cout<<"here1 "<<kpass<<endl;getchar();
//        cout<<"First: ";
//        for (int i=0;i<nume;i++)
//            cout<<eigs[i]<<" ";
//        cout<<endl;
        
        DiagonalizeS(stm,twosvec,eigs,kpass,nume);
//        cout<<"Second: ";
//        for (int i=0;i<nume;i++)
//            cout<<eigs[i]<<" ";
//        cout<<endl;
//        cout<<"here2"<<endl;getchar();
        
        *ierr=-abs(info);
        if (*ierr!=0) return;
        /* TeST for convergence on the absolute difference of eigenvalues between
         *  successive steps. Also SELect the unconverged eigenpairs and sort them
         *  by the largest magnitude in the last added NNCV rows of Svec. 	*/

        done=Test_Convergence(kpass,nume,neig,isel,twosvec,eigs,icv2,crite,critc,old,&nncv,incv2);
//        getchar();
//        cout<<"Done: "<<done<<endl;
        if ((done)||(kpass>=n)) {
        	omin_to=true;
            break;
        }
        if (kpass==lim) {
            /* Maximum size for expanding basis. Collapse basis, D, and S, Svec
             *  Consider the basis vectors found in TSTSEL for the newvec.	*/
            Collapse_Matrix(twosvec,b,ptb,lim);
            Collapse_Matrix(twosvec,ab,ptab,lim);
            Collapse_S(twosvec,stm,eigs,nume,lim);
            kpass=nume;
            
        }
//        for (int i=0;i<incv2.size();i++)
//            cout<<incv2[i]<<" ";
//        cout<<endl;
        /* Compute and add the new vectors. NNCV is set to the number of new
         *  vectors that have not converged. If none, DONE=true, exit.	*/
        Add_Vectors_To_Basis(nume,lim,mblock,kpass,critr,ortho,&nncv,incv2,diag->Diagonal(),twosvec,eigs,ab,b,icv2,&restart,&done);
        /*          ..An infinite loop is avoided since after a collapsing Svec=I
         *          ..=> Res=Di-lBi which is just computed and it is orthogonal.
         *          ..The following is to prevent an improbable infinite loop.  */
//        cout<<"Restart: "<<restart<<endl;
        if (!restart)
        {
            first=true;
        }
        else if (first)
        {
        	first=false;

            Collapse_Matrix(twosvec,b,ptb,kpass+nncv);
            Collapse_Matrix(twosvec,ab,ptab,kpass+nncv);
            Collapse_S(twosvec,stm,eigs,nume,lim);
            for (int i=0;i<nume;i++)
                eigs[i]=1.e+30;
            kpass=nume;
            continue;
        }
        else
        {
        	*ierr += 1024;
            omin_to=true;
            break;
        }
        if (done)
        {
        	omin_to=true;
            break;
        }
        //Add new columns in D and S, from the NNCV new vectors.
        Add_Vectors_To_AB_S(n,lim,hiend,kpass,nncv,b,ab,stm,diag);
        
        
        *nmv += nncv;
        kpass += nncv;
        ++nloops;
     
        if (userbreak) break;
        if (inf) maxiter=nloops+1;
    } while (nloops<maxiter || inf);
//    cout<<endl<<"It took: "<<nloops<<" ("<<maxiter<<" allowed max)"<<endl;
//    cout<<"Non converged: "<<nncv<<endl;
//    getchar();
    
    if (!omin_to)
    {
    	omin_to=false;
        kpass -= nncv;
    }

//    if (omin_to) {
//    	++nloops;
//    }

    /* Calculate final results. EIGVAL contains the eigenvalues, BASIS the
     *  eigenvectors, OLDVAL the eigenvalue differences, and SCRA1 residuals.	*/
    for (ntyp i=0; i<nume; i++)
        old[i]=fabs(old[i]-eigs[i]);
    
    if (kpass!=nume){
    Collapse_Matrix(twosvec,b,ptb,kpass);
    Collapse_Matrix(twosvec,ab,ptab,kpass);
    }
    /* i=1,NUME residual(i)= DCi-liBCi= newDi-linewBi                        ACPZ0579
     *  temporarily stored in AB(NUME*N+1)	*/
    
    for (ntyp i=0; i<nume; i++) {
        ab[nume]=ab[i];
        daxpy_vector(-eigs[i],b[i],ab[nume]);
        scra1[i]=sqrt(dot_vectors(ab[nume],ab[nume]));
    }
//    cout<<"edw"<<getchar();
    *loop=nloops;
    return;
}

void dscal(ntyp n, ftyp da, ftyp dx[], ntyp incx)
/*   scales a vector by a constant.
 uses unrolled loops for increment equal to one.
 jack dongarra, linpack, 3/11/78.
 modified 3/93 to return if incx .le. 0.
 modified 12/3/93, array(1) declarations changed to array(*) */
{
    if (n<=0 || incx<=0) return;
    if (incx!=1) {
    	//code for increment not equal to 1
    	ntyp nincx = n*incx;
        for (ntyp i=1; i<=nincx; i+=incx) dx[i -1] *= da;
        return;
    }
    /*      code for increment equal to 1
     
     clean-up loop  */
    ntyp m=n%5;
    if (m) {
    	for (ntyp i=1; i<=m; i++) dx[i -1] = da*dx[i -1];
        if (n<5) return;
    }
    for (ntyp i=m + 1; i<=n; i+=5) {
    	dx[i -1] *= da;
        dx[i] *= da;
        dx[i + 1] *= da;
        dx[i + 2] *= da;
        dx[i + 3] *= da;
    }
    return;
}

#endif


