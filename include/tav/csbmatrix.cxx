/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/**
 \author Alexander Volya  <http://www.volya.net>
 \author K. Kravvaris
 \file csbmatrix.cxx
 \brief Sparse matrix class, stores matrix in sparse blocks.
 */

#ifndef __Tav__CSBMatrix__
#define __Tav__CSBMatrix__

//io is not strictly neccessary, only when NDEBUG flag is not set.
#include <iostream>
#include <iomanip>
#include <cstdint>//for uint*_t
#include <fstream>//Read/Write from file
#include <cmath>//sqrt in ReadSparseMatrix
//omp.h is in seperate ifdef, the code should still compile without openmp.
#ifdef _OPENMP
#include <omp.h>
#endif

namespace tav
{
    
    class mutex_array
    {
    public:
        uint16_t n;//size of array
#ifdef _OPENMP
        omp_lock_t* mutexes;
#endif
        mutex_array(uint16_t n_):n(n_)
        {
#ifdef _OPENMP
            mutexes = new omp_lock_t[n];
#endif
        }
        ~mutex_array()
        {
#ifdef _OPENMP
            delete[] mutexes;
#endif
        }
        void init_all()
        {
#ifdef _OPENMP
            for (uint16_t i = 0; i < n; ++i)
                omp_init_lock(&(mutexes[i]));
#endif
        }
        void init(uint16_t i)
        {
#ifdef _OPENMP
            omp_init_lock(&(mutexes[i]));
#endif
        }
        void set(uint16_t i)
        {
#ifdef _OPENMP
            omp_set_lock(&(mutexes[i]));
#endif
        }
        bool test(uint16_t i)
        {
#ifdef _OPENMP
            return omp_test_lock(&(mutexes[i]));
#endif
            return true;
        }
        void unset(uint16_t i)
        {
#ifdef _OPENMP
            omp_unset_lock(&(mutexes[i]));
#endif
        }
    };
    //----------------End includes, start classes---------------------
    /*
     This class describes a single block in the CSB format.
     It contains:
     1 - isDiagonal: a Boolean variable to check whether this block belongs in the diagonal of
     the matrix (and multiply accordingly)
     2 - nnz: The number of non zero elements in this block.
     3 - indRow: An array storing the row index of all elements in the block. Indicies are
     Block specific not array specific.
     4 - indCol: An array storing the column index of all elements in the block. Indicies are
     Block specific not array specific.
     5 - elem: An array storing the elements of the block.
     
     Constructors:
     There are 3 constructors, a void one (to make pointers nicely), a bool one which just
     sets the diagonal and one that completely initializes the block by passing everything.
     
     Read/Write Functions:
     Read and Write functions from/to binary are provided.
     Note: these function take an already open stream as input and not a file name
     as they are only part of a total write function that occurs at the full
     matrix level, so it is much more helpful to pass this and continue writing
     wherever the stream has stopped. Obviously they can be used to write to
     seperate files too.
     
     WriteRow Functions:
     These are used to add elements to the block.
     WriteRowMove takes 2 arrays of indices (Row/Column) an array of elements and their size
     and appends them to the end of the old arrays.
     This is done by using std::move so old contents might not be there afterwards.
     
     WriteRowAt is similar however it writes 1 specific row, so does not need an array of
     row indices. This function exists to make reading files saved in other formats a bit faster.
     
     Block Multiplications:
     TimesVector multiplies the block with a vector
     TimesVectorTranspose multiplies the transpose of the block with a vector
     Note: no parallelism in multiplication exists for these 2 functions, nor is the
     vector initialized to 0 in any way as this is only part of the calculation.
     
     TODO: Read/Write from const char* filename
     Multiplications with Blocks of Vectors.
     
     */
    template<typename dtype = double, typename ubig_t = uint32_t, typename usmall_t = uint16_t>
    class CSBBlock
    {
    public:
        bool isDiagonal = false;//true if block lies along the diagonal.
        ubig_t nnz = 0;//number of non-zero elements in block
        usmall_t* indRow = nullptr;//row indices of all elements
        usmall_t* indCol = nullptr;//column indices of all elements
        dtype* elem = nullptr;//list of all elements.
        dtype* diagonal = nullptr;//stores the matrix elements that lie along the diagonal for fast access.
        
        CSBBlock():isDiagonal(false),nnz(0),indRow(nullptr),indCol(nullptr),elem(nullptr),diagonal(nullptr)
        {;}
        CSBBlock(bool isD):isDiagonal(isD),nnz(0),indRow(nullptr),indCol(nullptr),elem(nullptr),diagonal(nullptr){;}
        CSBBlock(ubig_t nnz_, usmall_t* indRow_,usmall_t* indCol_, dtype* elem_, bool isD = false):isDiagonal(isD),nnz(nnz_),indRow(nullptr),indCol(nullptr),elem(nullptr),diagonal(nullptr)
        {
            indRow = new usmall_t[nnz];
            indCol = new usmall_t[nnz];
            elem = new dtype[nnz];
            
            std::copy(indRow_,indRow_ + nnz, indRow);
            std::copy(indCol_,indCol_ + nnz, indCol);
            std::copy(elem_  ,elem_ + nnz  , elem  );
            
        }
        ~CSBBlock()
        {
            delete[] indRow;
            delete[] indCol;
            delete[] elem;
            delete[] diagonal;
        }
        usmall_t MakeDiagonal(usmall_t siz)
        {
            delete[] diagonal;
            diagonal = new dtype[siz];
            memset(diagonal, 0, siz * sizeof(dtype));
            usmall_t cnt = 0;
            for (ubig_t i = 0; i < nnz; ++i)
                if (indCol[i] == indRow[i])
                    diagonal[cnt++] = elem[i];
            return cnt-1;
        }
        
        //Write a simple extra row. assume it's a lower row for now.
        ubig_t WriteRowMove(ubig_t nnz_, usmall_t* indRow_,usmall_t* indCol_, dtype* elem_)
        {
            usmall_t *toDc = indCol, *toDr = indRow;
            dtype *toDe = elem;
            indCol = new usmall_t[nnz + nnz_];
            indRow = new usmall_t[nnz + nnz_];
            elem   = new dtype   [nnz + nnz_];
            std::move(toDr, toDr + nnz, indRow);
            std::move(indRow_,indRow_ + nnz_, indRow + nnz);
            std::move(toDc, toDc + nnz, indCol);
            std::move(indCol_,indCol_ + nnz_, indCol + nnz);
            std::move(toDe, toDe + nnz, elem);
            std::move(elem_,elem_ + nnz_, elem + nnz);
            nnz += nnz_;
            delete[] toDc;
            delete[] toDr;
            delete[] toDe;
            return nnz;
        }
        //Write a simple extra row. assume it's a lower row for now.
        ubig_t WriteRowAt(usmall_t rownum, ubig_t nnz_, usmall_t* indCol_, dtype* elem_)
        {
            usmall_t *toDc = indCol, *toDr = indRow;
            dtype *toDe = elem;
            indCol = new usmall_t[nnz + nnz_];
            indRow = new usmall_t[nnz + nnz_];
            elem   = new dtype   [nnz + nnz_];
            
            std::move(toDr, toDr + nnz, indRow);
            std::fill(indRow+nnz, indRow + nnz + nnz_ , rownum);
            
            std::move(toDc, toDc + nnz, indCol);
            std::move(indCol_,indCol_ + nnz_, indCol + nnz);
            std::move(toDe, toDe + nnz, elem);
            std::move(elem_,elem_ + nnz_, elem + nnz);
            nnz += nnz_;
            
            delete[] toDc;
            delete[] toDr;
            delete[] toDe;
            return nnz;
        }
        //assume block was empty. and write all. this is the same as constructor
        //moves data so it gets deleted.
        ubig_t WriteFullBlockMove(usmall_t* rows, usmall_t* cols, dtype* elems, ubig_t nnz_)
        {
            nnz = nnz_;
            indCol = new usmall_t[nnz];
            indRow = new usmall_t[nnz];
            elem = new dtype[nnz];
            std::move(cols, cols + nnz, indCol );
            std::move(rows, rows + nnz, indRow );
            std::move(elems  , elems   + nnz, elem);
            
            
            return nnz;
        }
        int Read(std::ifstream& inpf)
        {
            inpf.read((char*)(&isDiagonal),sizeof(bool)); //we have saved n
            inpf.read((char*)(&nnz)  ,sizeof(ubig_t));
            
            indCol = new usmall_t[nnz];
            indRow = new usmall_t[nnz];
            elem = new dtype[nnz];
            
            inpf.read((char*)(indCol),nnz*sizeof(usmall_t));
            inpf.read((char*)(indRow),nnz*sizeof(usmall_t));
            inpf.read((char*)(elem)  ,nnz*sizeof(dtype));
            
            return 0;
        }
        
        int Write(std::ofstream& outf)
        {
            outf.write((char*)(&isDiagonal),sizeof(bool)); //we have saved n
            outf.write((char*)(&nnz)  ,sizeof(ubig_t));
            outf.write((char*)(indCol),nnz*sizeof(usmall_t));
            outf.write((char*)(indRow),nnz*sizeof(usmall_t));
            outf.write((char*)(elem)  ,nnz*sizeof(dtype));
            return 0;
        }
        
        template< class cnumber>
        int TimesVector(cnumber* y, cnumber* x)
        {
            if (isDiagonal)
            {
                for (ubig_t i = 0; i < nnz; ++i)
                {
                    y[indRow[i]] += elem[i] * x[indCol[i]];
                    if (indRow[i] != indCol[i])
                        y[indCol[i]] += elem[i] * x[indRow[i]];
                }
                
            }
            else
            {
                for (ubig_t i = 0; i < nnz; ++i)
                {
                    y[indRow[i]] += elem[i] * x[indCol[i]];
                }
            }
            return 0;
        }
        
        template< class cnumber>
        int TimesMatrix(cnumber** y, cnumber** x, usmall_t m)
        {
            if (isDiagonal)
            {
                for (usmall_t im = 0; im < m; ++im)
                {
                    for (ubig_t i = 0; i < nnz; ++i)
                    {
                        y[im][indRow[i]] += elem[i] * x[im][indCol[i]];
                        if (indRow[i] != indCol[i])
                            y[im][indCol[i]] += elem[i] * x[im][indRow[i]];
                    }
                }
                
            }
            else
            {
                for (usmall_t im = 0; im < m; ++im)
                {
                    for (ubig_t i = 0; i < nnz; ++i)
                    {
                        y[im][indRow[i]] += elem[i] * x[im][indCol[i]];
                    }
                }
            }
            return 0;
        }
        
        template< class cnumber>
        int TimesVectorTranspose(cnumber* y, cnumber* x)
        {
            for (ubig_t i = 0; i < nnz; ++i)
            {
                y[indCol[i]] += elem[i] * x[indRow[i]];
            }
            return 0;
        }
        template< class cnumber>
        int TimesMatrixTranspose(cnumber** y, cnumber** x, usmall_t m)
        {
            for (usmall_t im = 0; im < m; ++im)
            {
                for (ubig_t i = 0; i < nnz; ++i)
                {
                    y[im][indCol[i]] += elem[i] * x[im][indRow[i]];
                }
            }
            return 0;
        }
#ifndef NDEBUG
        double sizeMB()
        {
            size_t size_in_bytes = sizeof(bool) + sizeof(ubig_t) + nnz * (2*sizeof(usmall_t) + sizeof(dtype));
            return size_in_bytes * 1.e-6;
        }
#endif
    };
    
    
    template<typename dtype = double, typename ubig_t = uint32_t, typename usmall_t = uint16_t>
    class SparseMatrixCSB
    {
    public:
        bool isSymmetric;//true if matrix is symmetric
        ubig_t n;//size of matrix
        usmall_t BlockRow;//number of rows in a block
        usmall_t BlockCol;//number of columns in a block
        usmall_t numBlockRow;//number of rows of blocks
        usmall_t numBlockCol;//number of columns of blocks
        
        CSBBlock<dtype, ubig_t, usmall_t>** Blocks;//The actual blocks.
        
        SparseMatrixCSB():isSymmetric(false),n(0),
        BlockRow(0),BlockCol(0),numBlockRow(0),numBlockCol(0),Blocks(nullptr)
        {;}
        
        SparseMatrixCSB(usmall_t BR, usmall_t BC, bool isS=false):isSymmetric(isS),n(0),
        BlockRow(BR),BlockCol(BC),numBlockRow(0),numBlockCol(0),Blocks(nullptr)
        {;}
        SparseMatrixCSB(ubig_t n_,usmall_t BR, usmall_t BC, bool isS=false):isSymmetric(isS),n(n_),
        BlockRow(BR),BlockCol(BC),numBlockRow(0),numBlockCol(0),Blocks(nullptr)
        {
            numBlockRow = n/BlockRow;
            numBlockCol = n/BlockCol;
            if (numBlockRow * BlockRow != n) numBlockRow++;
            if (numBlockCol * BlockCol != n) numBlockCol++;
            
            Blocks = new CSBBlock<dtype, ubig_t, usmall_t>* [numBlockRow];
            for (usmall_t i = 0; i < numBlockRow; ++i)
            {
                Blocks[i] = new CSBBlock<dtype, ubig_t, usmall_t> [numBlockCol];
                if (isSymmetric) Blocks[i][i].isDiagonal = true;
                for (usmall_t j = 0; j < numBlockCol; ++j)
                    Blocks[i][j].nnz = 0;
            }
            

        }
        ~SparseMatrixCSB()
        {
            for (ubig_t i = 0; i < numBlockRow ; ++i)
                delete[] Blocks[i];
            delete[] Blocks;
        }
        usmall_t MakeDiagonals()
        {
            for (usmall_t i = 0; i < numBlockRow; ++i)
            {
                if (Blocks[i][i].isDiagonal)
                    Blocks[i][i].MakeDiagonal(BlockRow);//assume blocks are square.
            }
            return 0;
        }
        ubig_t ReadDiagonal(dtype *x)
        {
            
//            for (ubig_t i = 0; i < numBlockCol ; ++i)
//            {
//                //            std::copy(x + i * BlockRow,x + (i + 1) *BlockRow, Blocks[i][i].diagonal);
//                for (usmall_t j = 0; j < BlockRow; ++j)
//                    x[i * BlockRow + j] = Blocks[i][i].diagonal[j];
//            }
            memset(x, 0, n*sizeof(dtype));
//            return 0;
            for (ubig_t i = 0; i < this->size(); ++i)
                x[i] = Blocks[i/BlockRow][i/BlockRow].diagonal[i % BlockRow];
            
//            for (usmall_t i = 0; i < numBlockRow; ++i)
//            {
//                for (ubig_t j = 0; j < Blocks[i][i].nnz; ++j)
//                {
//                    if (Blocks[i][i].indRow[j] == Blocks[i][i].indCol[j])
//                    {
//                        std::cout << "writing: " << i*BlockRow << " " <<  Blocks[i][i].indRow[j] << " " << Blocks[i][i].indRow[j] << std::endl;
//                        x[i*BlockRow + Blocks[i][i].indRow[j]] = Blocks[i][i].elem[j];
//                    }
//                }
//            }
            return n;
        }
        
        ubig_t WriteRowMove(ubig_t rownum,dtype* elem, ubig_t* ColInds,ubig_t nz)
        {
            //we have a full row, we need to break it apart in smaller ones
            //to write to every block.
            ubig_t nnzi = 0;
            usmall_t blRow = rownum/BlockRow;
            
            ubig_t* cls = new ubig_t[numBlockCol];
            usmall_t* colInd = new usmall_t[nz];
            memset(cls,0,numBlockCol * sizeof(ubig_t));//initialize to 0
            
            for (ubig_t i = 0; i < nz; ++i)
            {
                cls[ColInds[i] / BlockCol]++;
                colInd[i] = usmall_t(ColInds[i] % BlockCol);
            }
            for (usmall_t i = 0; i < numBlockCol; ++i)
            {
                Blocks[blRow][i].WriteRowAt(usmall_t(rownum % BlockRow), cls[i], colInd + nnzi, elem+nnzi);
                nnzi += cls[i];
            }
            
            delete[] cls;
            delete[] colInd;
            return nz - nnzi;
            
        }
        int ReadFromSparseFile(const char* filename)
        {
            std::ifstream outf;
            outf.open(filename, std::ios::binary);
            
            if (!outf)   {               // if the file does not exist
                std::cerr<<"File "<<filename<<" is not found"<<std::endl; // return error message
                return 0;
            }
            //TODO: Delete old stuff if any
            std::cout << "Reading from file: " << filename << std::endl;
            
            typedef unsigned long long uint_nbs;
            uint_nbs n_ = 0;//THIS IS HERE FOR COMPATIBILITY WITH SPARSE MATRIX.
            //Read now new things
            outf.read((char*)(&n_),sizeof(uint_nbs)); //got n
            
            n = ubig_t(n_);
            
            std::cout << "Read Sparse file with n = " << n << std::endl;
            BlockRow = std::sqrt(double(n)) + 1;
            BlockCol = BlockRow;
            
            numBlockRow = n/BlockRow;
            numBlockCol = n/BlockCol;
            if (numBlockRow * BlockRow != n) numBlockRow++;
            if (numBlockCol * BlockCol != n) numBlockCol++;
            
            Blocks = new CSBBlock<dtype, ubig_t, usmall_t>* [numBlockRow];
            for (usmall_t i = 0; i < numBlockRow; ++i)
            {
                Blocks[i] = new CSBBlock<dtype, ubig_t, usmall_t> [numBlockCol];
                for (usmall_t j = 0; j < numBlockCol; ++j)
                    Blocks[i][j].nnz = 0;
            }
            
            //make them diagonal.
            for (usmall_t i = 0; i < numBlockRow; ++i) Blocks[i][i].isDiagonal = true;
            
            
            unsigned nz;
            uint_nbs* ind;
            ubig_t* inbig;
            dtype* elem;
            
            for(ubig_t i=0; i < n; i++) {
                if (i%100 == 0) std::cout << "\rDone Reading Sparse: " << double(i)/double(n) * 100. << "%     " << std::flush;
                outf.read((char*)(&(nz)),sizeof(unsigned));//read non-zero stuff
                ind  = new uint_nbs [nz];
                elem = new dtype [nz];
                inbig= new ubig_t [nz];
                outf.read((char*)(ind),nz*sizeof(uint_nbs)); //got ind
                outf.read((char*)(elem),nz*sizeof(dtype)); //got elem
                
                for (ubig_t j = 0; j < nz; ++j)
                    inbig[j] = ubig_t(ind[j]);
                //
                this->WriteRowMove(i, elem, inbig, ubig_t(nz));
                delete[] ind;
                delete[] elem;
                delete[] inbig;
                //
            }
            std::cout << std::endl;
            outf.close();
            MakeDiagonals();
            return 0;
            
        }
        int Read(const char* file)
        {

            //open stream
            std::ifstream inpf;
            inpf.open(file,std::ios::binary);
            if (!inpf)   {               // if the file does not exist
                std::cerr<<"File "<<file<<" is not found"<<std::endl; // return error message
                return 0;
            }
            //TODO: Delete old stuff if any
            std::cout << "Reading from file: " << file << std::endl;

            
            inpf.read((char*)(&isSymmetric),sizeof(bool)); //we have read the bool
            inpf.read((char*)(&n),sizeof(ubig_t)); //we have read n
            //read block related sizes
            inpf.read((char*)(&BlockRow),sizeof(usmall_t));
            inpf.read((char*)(&BlockCol),sizeof(usmall_t));
            inpf.read((char*)(&numBlockRow),sizeof(usmall_t));
            inpf.read((char*)(&numBlockCol),sizeof(usmall_t));
            
            //make empty block array of correct size
            Blocks = new CSBBlock<dtype, ubig_t, usmall_t>* [numBlockRow];
            for (usmall_t i = 0; i < numBlockRow; ++i)
                Blocks[i] = new CSBBlock<dtype, ubig_t, usmall_t> [numBlockCol];
            
            //read each element seperately
            for (usmall_t i = 0; i < numBlockRow; ++i)
            {
                for (usmall_t j = 0; j < numBlockCol; ++j)
                    Blocks[i][j].Read(inpf);
            }
            MakeDiagonals();
            return 0;
        }
        int Write(const char* file)
        {
            //open Stream
            std::ofstream outf;
            outf.open(file, std::ios::binary);
            
            outf.write((char*)(&isSymmetric),sizeof(bool)); //we have saved the symmetricity
            outf.write((char*)(&n),sizeof(ubig_t)); //we have saved n
            //save all block size stuff
            outf.write((char*)(&BlockRow),sizeof(usmall_t));
            outf.write((char*)(&BlockCol),sizeof(usmall_t));
            outf.write((char*)(&numBlockRow),sizeof(usmall_t));
            outf.write((char*)(&numBlockCol),sizeof(usmall_t));
            
            //save each block seperately
            for (usmall_t i = 0; i < numBlockRow; ++i)
            {
                for (usmall_t j = 0; j < numBlockCol; ++j)
                    Blocks[i][j].Write(outf);
            }
            
            
            return 0;
        }
        template<class cnumber>
        int TimesVector(cnumber* y, cnumber* x)
        {
            memset(y,0,n*sizeof(cnumber));
            if (isSymmetric)
            {
                mutex_array mutexes(numBlockCol);
                mutexes.init_all();
#ifdef _OPENMP
#pragma omp parallel for
#endif
                for (usmall_t j = 0; j < numBlockRow; ++j)
                {
                    while(!mutexes.test(j));
                    for (usmall_t i = 0; i < numBlockCol; ++i)
                    {
                        
                        if (i < j)
                            Blocks[i][j].TimesVectorTranspose(y + j * BlockRow, x + i * BlockCol);
                        else
                            Blocks[j][i].TimesVector(y + j * BlockRow, x + i * BlockCol);
                        
                    }
                    mutexes.unset(j);
                }
            }
            else
            {
                for (usmall_t j = 0; j < numBlockRow; ++j)
                    for (usmall_t i = 0; i < numBlockCol; ++i)
                        Blocks[j][i].TimesVector(y + j * BlockRow, x + i * BlockCol);
            }
            
            return 0;
        }
        int TimesMatrix2(dtype** y, dtype** x, usmall_t m)
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (usmall_t im = 0; im < m; ++im)
                this->TimesVector(y[im],x[im]);
            return 0;
        }
        //ok the idea here is we pre block the array of input- output vectors
        //so that blocks hopefully stay in cache longer.
        int TimesMatrix(dtype** y, dtype** x, usmall_t m) const
        {
            dtype*** xpoint = new dtype** [numBlockCol];
            dtype*** ypoint = new dtype** [numBlockCol];
            for (usmall_t i = 0; i < m; ++i) memset(y[i],0,n*sizeof(dtype));
            for (usmall_t j = 0; j < numBlockCol; ++j)
            {
                ypoint[j] = new dtype*[m];
                xpoint[j] = new dtype*[m];
                for (usmall_t i = 0; i < m; ++i)
                {
                    xpoint[j][i] = &(x[i][j*BlockCol]);
                    ypoint[j][i] = &(y[i][j*BlockCol]);
                }
                
            }
            if (isSymmetric)
            {
                mutex_array mutexes(numBlockCol);
                mutexes.init_all();
#ifdef _OPENMP
#pragma omp parallel for
#endif
                for (usmall_t j = 0; j < numBlockRow; ++j)
                {
                    while(!mutexes.test(j));
                    for (usmall_t i = 0; i < numBlockCol; ++i)
                    {
                        
                        if (i < j)
                            Blocks[i][j].TimesMatrixTranspose(ypoint[j], xpoint[i],m);
                        else
                            Blocks[j][i].TimesMatrix(ypoint[j], xpoint[i],m);
                        
                    }
                    mutexes.unset(j);
                }
            }
            
            return 0;
        }
        int TimesMatrix(dtype* y, dtype* x, usmall_t m) const
        {
            dtype **yy=new dtype*[m]; //pointers to vectors
            dtype **xx=new dtype*[m]; //pointers to out vectors
            for(unsigned im=0;im<m;im++) {//loop over m vectors
                yy[im]=y+n*im; //shift pointer so yy[i] is a column vector
                xx[im]=x+n*im;}
            TimesMatrix(yy, xx, m);
            delete [] yy;
            delete [] xx;
        return 0;
        }
        int TimesMatrix3(dtype** y, dtype** x, usmall_t m)
        {
            dtype*** xpoint = new dtype** [numBlockCol];
            dtype*** ypoint = new dtype** [numBlockCol];
            for (usmall_t i = 0; i < m; ++i) memset(y[i],0,n*sizeof(dtype));
            for (usmall_t j = 0; j < numBlockCol; ++j)
            {
                ypoint[j] = new dtype*[m];
                xpoint[j] = new dtype*[m];
                for (usmall_t i = 0; i < m; ++i)
                {
                    xpoint[j][i] = &(x[i][j*BlockCol]);
                    ypoint[j][i] = &(y[i][j*BlockCol]);
                }
                
            }
            if (isSymmetric)
            {
                mutex_array mutexes(numBlockCol);
                mutexes.init_all();
#ifdef _OPENMP
#pragma omp parallel for
#endif
                for (usmall_t j = 0; j < numBlockRow; ++j)
                {
                    
                    while(!mutexes.test(j));
                    for (usmall_t i = j; i < numBlockCol; ++i)
                    {
                        if (i != j)
                            while(!mutexes.test(i));
                        Blocks[i][j].TimesMatrix(ypoint[j],xpoint[i],m);
                        Blocks[i][j].TimesMatrixTranspose(ypoint[i],xpoint[j],m);
                        if (i!=j)
                            mutexes.unset(i);
                    }
                    mutexes.unset(j);
                }
            }
            
            return 0;
        }
        
        size_t size()
        {
            return size_t(n);
        }
        //these are here just for debugging
        //they print the number of non-zero elements
        //per block. Nothing to see here, move along.
#ifndef NDEBUG
        void PrintZeroStructure()
        {
            for (usmall_t i = 0; i < numBlockRow; ++i)
            {
                for (usmall_t j = 0; j < numBlockCol; ++j)
                    std::cout << Blocks[i][j].nnz << " ";
                std::cout << std::endl;
            }
        }
        void PrintZeroStructure2()
        {
            int nzq = 0;
            for (usmall_t i = 0; i < numBlockRow; ++i)
            {
                nzq = 0;
                for (usmall_t j = 0; j < numBlockCol; ++j)
                    nzq += Blocks[i][j].nnz;
                std::cout << nzq <<  std::endl;
            }
        }
        double sizeMB()
        {
            size_t correction = 4*sizeof(usmall_t) + sizeof(bool) + sizeof(ubig_t);
            double sizeM = 0.;
            for (usmall_t i = 0; i < numBlockRow; ++i)
            {
                for (usmall_t j = 0; j < numBlockCol; ++j)
                    sizeM += Blocks[i][j].sizeMB();
            }
            return sizeM + correction * 1.e-6;
        }
#endif
    };
    
};//end namespace tav
#endif 
