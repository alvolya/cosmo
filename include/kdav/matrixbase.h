
#ifndef __MATRIX_BASE
#define __MATRIX_BASE


namespace MatB {
    
class MatrixBase
{
	public:
        // the method returns pointer for the vector of the diagonal elements
    	virtual ftyp *Diagonal() const=0;

        /* computes the product C of matrix A (of size n) and a block m of
           vectors B (of size n), that is C=A B*/
        virtual void TimesMatrix(ftyp *, ftyp *, ntyp ) const=0;
        virtual void TimesMatrix(ftyp **, ftyp **, ntyp ) const=0;

        /* returns the dimension of the table containing the diagonal elements
           (equal to dimension of matrix A) */
        virtual ntyp Size() const=0;

        //this operator returns  the value of the i-th diagonal element
        ftyp operator[](const ntyp &i) const {return *(Diagonal()+i);};
};
}
#endif

