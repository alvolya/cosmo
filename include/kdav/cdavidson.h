
#ifndef _CDAVIDSON_H
#define _CDAVIDSON_H

//#include "gridparam.h"
#include "defaulttd.h"
#include "matrixbase.h"
#include "fsdavidson.cxx"
// Structure stores eigenvalue, its residual and
// difference beetwen two final iterations
struct eig_val{
	ftyp Eigenvalue;
    ftyp Eigval_differences;
    ftyp Residual;
};

class Davidson 
{
    //private : i needed n and work to be public
    public :
    	//number of eigenpair to be found
        ntyp neig,
        
             /* An integer denoting the completions status:
		   		if IERR = 1024 then orthogonalization failured
		   		if IERR = 2048 then the maximum number of iterations
           						has benn reached
           		if IERR = 4096 then the algorithm has been broken by user
           					(press ENTER and y)
           		Becouse this class uses service exceptions mechanism,
           		the other errors are submitted by means of this mechanism */
        	 ierr,

             /* The index of the lowest eigepair to be computed. If
            	(ILOW<=0).or.(ILOW>N), the selected eigenpairs
            	to be computed should be contained in array ISELEC.
            	(Modified on exit). */
        	 ilow,

            /*	Number of Initial Vector estimates provided by the user.
            	If NIV is in the range:  (NUME).LE.(NIV).LE.(LIM),
            	the first NIV columns of size N of WORK should contain
            	the estimates (see below). In all other cases of NIV,
            	the program generates initial estimates.
            	This implementation of Davidson method does not permit
				to use Initial Vector estimates provided by the user. */
             niv,

             /*	The number of Matrix-vector(M-V) multiplies. Each matrix
            	reference can have up to size(block) M-V multiplies.  */
             nmv,

             /*	The index of the highest eigenpair to be computed.
            	Considered ONLY when ILOW is in the range
            	(0<ILOW<=N). (Modified on exit). */
             ihigh;

        	 /*	The size of the integer workspace. It must be at least
            	as large as:
                                 6*LIM + NUME

            	If LIM or NUME needs to be increased, the space should
            	also be increased accordingly. For given IWRSZ and
            	IIWSZ one can calculate how big a problem one can
            	solve (LIM,NUME).   */
        ntyp iiwsz,

             /* The upper limit on the dimension of the expanding basis.
             	NUME.LT.LIM.LE.N must hold. The case LIM=NUME is allowed
            	only for LIM=NUME=N. The choice of LIM depends on the
            	available workspace (see below). If the space is
            	available it is preferable to have a large LIM, but not
            	larger than NUME$+$40.
            	In this implementation of Davidson method LIM=LIMMAX */
        	 lim,

             /*	The size of the real workspace. It must be at least as
            	large as:
                    2*N*LIM + LIM*LIM + (NUME+10)*LIM + NUME */
             irwsz;

		/* An object diag is a MatrixBase class. The MatrixBase class is
		   responsible for diagonal values matrix (calculates
           it and makes it available with use of [] operator)
           and provides method computes the product of matrix
           A and a block of vectors B(n,m), C=A B where A(nxn)  */
        const MatB::MatrixBase *diag;

        char initfromfiles[255];
        ntyp nume,
        	 n, // Order of the matrix A  (Ax=Ex).
             numemax, limmax;

             /* Array of size LIM holding the user specified indices
            	for the eigenpairs to be computed. Considered only when
            	(ILOW<=0)||(ILOW>N). The indices are read from
            	the first position until a non positive integer is met.
               	Example: if N=500, ILOW=0, and ISELEC(1)=495,
               	ISELEC(2)=497, ISELEC(3)=-1, the program will find
               	2 of the highest eigenpairs, pairs 495 and 497.
            	Any order of indices is acceptable (Modified on exit). */
    ntyp *iselec;

             /*	ntyp work array of size IIWSZ. Used as scrath array
            	for indices and for use in the LAPACK routines.  */
        	/* *iwork*/

        	/* WORK(1)
        			The first NUME*N locations contain the approximations
        			to the NUME extreme eigenvectors. If the lowest eigenpairs
            		are required, (HIEND=false), eigenvectors appear in
            		ascending order, otherwise (HIEND=false), they appear in
            		descending order. If only some are requested, the order
            		is the above one for all the NUME extreme eigenvectors,
            		but convergence has been reached only for the selected
            		ones. The rest are the current approximations to the
            		non-selected eigenvectors.

        		WORK(NUME*N+1)
            		The next NUME locations contain the approximations to
            		the NUME extreme eigenvalues, corresponding to the above
            		NUME eigenvectors. The same ordering and convergence
            		status applies here as well.

				WORK(NUME*N+NUME+1)
            		The next NUME locations contain the corresponding values
            		of ABS(EIGVAL-VALOLD) of the NUME above eigenvalues, of
            		the last step of the algorithm.

				WORK(NUME*N+NUME+NUME+1)
            		The next NUME locations contain the corresponding
            		residual norms of the NUME above eigenvectors, of the
            		last step.  */
    vector<vector<double> > twosvec,b,ab;
    vector<double> smt;
    vector<ftyp> eigs,scratch1,old;
    vector<ntyp> icv,scra2,scra3;

             /* Logical. If true on exit the highest eigenpairs are
            	found in descending order. Otherwise, the lowest
            	eigenpairs are arranged in ascending order.  */
        bool hiend,

        	 tab, Init,Init2;

        void InitVectors();
        /*Sets default values following field:
            	ierr=0;
    			crite = 1e-15; critc = 1e-12; critr = 1e-8; ortho = 1e-9;
                niv = 0;
    			mblock = 1;
    			numemax = Nume;
    			maxiter = MAX( numemax*40, 200 );
    			nmv=0;
    			neig = 0;
    			limmax=numemax+20;
    			iiwsz = 6*limmax + numemax;
    			lim = limmax;
    			irwsz = 2*n*limmax + limmax*limmax + (numemax+10)*limmax + numemax;
    			hiend=false;
    			norm=0;
    			silence=false;
    			loop=0;  */
        void SetDefault(ntyp);

        	/* Check parameters setting in order to detect errors.
            	This method submittes exceptions ErrDavidson class. */
        void CheckError(ntyp iselec[]);

//        ftyp suma(ntyp n, ftyp tab[]);
        	 //Returns sum of elements of tab vector dimension n.

        	/* Normalizetes tab vector dimension n to:
            	1. in case of 1D : dx
                2. in case of 2D : dx*dy
                3. in case of 3D : dx*dy*dz */
//        ntyp normalizacja(ntyp n, ftyp tab[]);

        	//Execute Davidson method.
        void ExecuteDavidson();

   // public :
             /* Convergence threshold for eigenvalues.
            	If ABS(EIGVAL-VALOLD) is less than CRITE for all wanted
            	eigenvalues, convergence is signaled. */
        ftyp crite,

             /* Convergence threshold for the coefficients of the last
            	added basis vector(s). If all of those corresponding to
            	unconverged eigenpairs are less than CRITC convergence
            	is signaled. */
        	 critc,

             /*	Convergence threshold for residual vector norms. If
            	all the residual norms ||Ax_i-l_ix_i|| of the targeted
            	x_i are less than CRITR convergence is signaled.
            	If ANY of the criteria are satisfied the algorithm stops */
             critr,

             /*	The threshold over which loss of orthogonality is
            	assumed. Usually ORTHO.LE.CRITR*10 but the process can
            	be skipped by setting ORTHO to a large number(eg,1.D+3). */
             ortho,

             sum;

        	 /* Upper bound on the number of iterations of the
            	algorithm. When MAXITER is exceeded the algorithm stops.
            	A typical MAXITER can be MAX(200,NUME*40), but it can
            	be increased as needed. */
        ntyp maxiter,

             /*	Number of vectors to be targeted in each iteration.
            	1<=MBLOCK<=(No. EiGenpairs wanted) should hold.
            	Large block size reduces the number of iterations
            	(matrix acceses) but increases the matrix-vector
            	multiplies. It should be used when the matrix accese
            	is expensive (disc, recomputed or distributed). */
        	 mblock,

             //if vectors are normalized then norm=1 else norm=0
             norm,
             /* The number of iterations it took to reach convergence.
           		This is also the number of matrix references. */
             loop;

             /*	if silence=0 an additional information is displayed */
        bool silence,

             /* If true then upper bound on the numeber of iterations
                is equal to maximum number of type ntyp */
        	 inf;

        /* 	The eigvalmax parameter is the maximum number of
            searche  eigenvalue and eigvalmin is the minimum
            number of searched eigenvalue. The program finds
            eigenvalues from range [eigvalmin, eigvalmax].
            If eigvalmin=0 the program produces only eigenvalue of
            number eigenvalmax (for eigenvalues < eigenvalmax
            convergence has not been reached). */
        Davidson(const MatB::MatrixBase *, ntyp eigvalmax=1, ntyp eigvalmin=0);

        /* 	In this constructor the second argument is a table
            of indices for the eigenpairs to be computed.
            Iselect must have 0 (or number<0) for last positions.
            For example, if iselect={1,5,7,0} program finds
            eigenpairs of indices 1, 5 and 7. */
        Davidson(const MatB::MatrixBase *, ntyp *);

        // Starts the calculations.
        void Calculate();

        /* reads the initial vectors from files. The first argument is a file
           name, which includes paths name to text files containing the initial
           vectors  */
		void ReadInitVectors(const char *);

        // This method sets the maximum index of eigenpairs to be found.
        void SetNumeMax(ntyp );

		/*	This method sets the upper limit on the order of
            the expanding basis. It must be used after SetNumeMax
            (in case one uses SetNumeMax), which restores the default
            limit on size of the expanding basis. */
        void SetLimMax(ntyp );

        /* The argument is a table of indices for the eigenpairs to be computed.
           The last value in this table must be equal or less than zero, e.g.
           {1,5,7,0} (eigenpairs of indices 1, 5 and 7)*/
        void SetEigSelect(ntyp *);

        //sets the Number of Initial Vector estimates provided by the user
        void SetNiv(ntyp );

        //gets the maximum index of the eigenpair to be found
        ntyp GetNumeMax() const;

        //gets the upper limit on the order of the expanding basis
        ntyp GetLimMax() const;

        /* this method returns 2048 (if the numer of executed iterations has
        exceeded the maximum number of iterations), 4096 (if the program has
        been stopped by the user) or 1024 (if orthogonalization failed) */
        ntyp GetIerr() const; //Returns ierr value (1024, 2048 or 4096)

        bool GetHiend() const; //Returns hiend value
        ntyp GetMatrixVectProd() const; //Returns nmv value

        /* 	
            This methods return residual of n-th eigenvector and save the
            components of an eigenvector to the memory block pointed by
            parameter vector.  If the eigenvector needs to be normalized,
            norm=1 should be set, else norm=0 (normalization is made to grid
            step and not to 1, this means that the sum of all vector elements
            gives 1). If normalization has failed then the value of the field
            norm of Davidson class is set to zero. */
//			 ftyp EigenFunction(ntyp n, ftyp   *, ntyp);
    

        /*	saves to structure of type eig_val the following information about
        the n-th eigenvalue: its value, the difference of eigenvalues from the
        last two iterations and the residual    */
        void EigenValue(eig_val &, 
		ntyp n ///<eigenvalue number counting from zero
		) const;
		void Save(std::ofstream &outputf, 
		ntyp nx 
		) const;
		void Save(std::ofstream &outputf) const;
        void Save(const char*) const; //<save all vectors to file
		void Save(const char*, ntyp n) const; //<save some number of vectors to file
		void SaveVector(std::ofstream &outputf, ntyp) const; //<save individual pair-vector
		void PrintEigenValue(ntyp nx) const;
		void PrintEigenValue() const;
        ~Davidson();

};

#endif //_CDAVIDSON_H

