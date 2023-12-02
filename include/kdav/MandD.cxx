/*! \file MandD.cxx
\ingroup gp_nr
*/
#ifndef __MANDD__CXX__
#define __MANDD__CXX__
//helper class that only adds diagonal to a matrix which is stored
template<typename SomeMatrix>
class MyMatrixDiagonal: public MatB::MatrixBase, SomeMatrix
{
        public:
		int dim;
		ftyp *dd;
		MyMatrixDiagonal(const SomeMatrix &MX) : SomeMatrix(MX) {
		dim=SomeMatrix::size();
		dd=new ftyp [dim];
		SomeMatrix::ReadDiagonal(dd);
		}
		~MyMatrixDiagonal() {delete [] dd;}  
        // the method returns pointer for the vector of the diagonal elements
        ftyp *Diagonal() const {return dd;};

        /* computes the product C of matrix A (of size n) and a block m of
           vectors B (of size n), that is C=A B*/

		   
        void TimesMatrix(ftyp *y, ftyp *x, ntyp m)  const
		{ 
		SomeMatrix::TimesMatrix(y,x,int(m));
		};
        void TimesMatrix(ftyp **y, ftyp **x, ntyp m)  const
        {
            SomeMatrix::TimesMatrix(y,x,int(m));
        };
		

        /* returns the dimension of the table containing the diagonal elements
           (equal to dimension of matrix A) */
        ntyp Size() const {return dim;};

        //this operator returns  the value of the i-th diagonal element
        ftyp operator[](const ntyp &i) const {return *(Diagonal()+i);};
};

#endif //__MANDD__CXX__
