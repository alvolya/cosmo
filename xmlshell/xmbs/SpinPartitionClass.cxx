/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


#ifndef __SPIN__PARTITION__CLASS__CXX__
#define __SPIN__PARTITION__CLASS__CXX__

#include <vector>
#include <iostream>
using namespace std;
//#include <tav/sort.cxx> //Locate is not needed, see below
/** This is a class that spin permutations
We always use projections so that \f$ \tilde{m}=m+j \f$
\todo Class is more general then just partition list, if dublicates are
allowed then construction of m-partitions is not needed.
*/
class SpinPartition: public vector<int *> {
private:
    int *partition;
    int w; ///<asset size
    int *asset; ///< ordered assets to be summed
    void Initiate(int M,int N, int W)
    {
        parts=N;
        sum=M; //shift the value of the sum, each particle gets J
        w=W;
        if (parts!=0) partition=new int [parts];
        asset=new int [w];
        for (int i=0;i<w;i++) {asset[i]=i;}
    }
    bool Make_(
               int parts_, ///<number of parts
               int sum_,  ///< that above partitions add upto
               int w_ ///<asset size
    )
    {
        //static int initparts=parts;
        if (parts_==0)
        {
            if (sum_!=0)
                return false;//
            else
            {
                int *tmp;
                push_back(tmp);
                return true;
            }
        }
        //int i=tav::Locate(w_,asset,sum_); ///try to allocate entire sum to one element
        int i=((sum_ >= w_) ? w_ - 1 : sum_); //this works because asset=0,1,2,3 array and locating there is trivial sum_>=0
        if (parts_==1)
        {
            if ((i==-1)||(asset[i]!=sum_))
                return false;
            else
            { //successful partition
                partition[parts_-1]=i;
//                {
                    int *tmp =new int [parts];
                    for (int q=0;q<parts;q++)
                    {
                        tmp[q]=partition[q];
//                        cout<<partition[q]<<" ";
                    }
//                cout<<endl;
                
                    push_back(tmp);
//                     delete[] tmp; //(?) this wasn't here before. tmp never deleted? why braces?
//                }
                return true;
            }
        } ///end if 1
        
        for (int j=0; j<=i; j++) {
            ///assume last number is j
            partition[parts_-1]=j;
            Make_(parts_-1,sum_-asset[j],j);
        }
        return false;
    }
    
public:
    int sum;
    int parts;
    int J;
    SpinPartition(int M,int N, int JJ) {
        Initiate(M,N,JJ);
    }
    
    SpinPartition() {
        sum=0;
        parts=0;
    }
    
    ~SpinPartition(void){
        delete [] asset;
        if (parts!=0) {
            delete [] partition;
            iterator ii;
            for (ii=begin(); ii!=end();ii++)  delete [] *ii;
        }
    }
    int Make(int M,int N, int JJ)
    {
        Initiate(M,N,JJ);
        Make_(parts,sum,w);
        return size();
    }
    int Make(void)
    {
        Make_(parts,sum,w);
        return size();
    }
    
    friend std::ostream& operator << (std::ostream &stream,
                                      SpinPartition &P) {
        for (SpinPartition::iterator ii=P.begin(); ii!=P.end();ii++)
        {for (int i=0;i<P.parts;i++) 
            stream<<(*ii)[i]<<" ";
            stream<<'\n';
        }
        return stream;
    }
    //   int Make(void) {Make_(sum,parts); return size();}
    /*   
     friend std::ostream& operator << (std::ostream &stream,
     Partition &P) { 
     for (Partition::iterator ii=P.begin(); ii!=P.end();ii++) 
     {for (int i=0;i<P.parts;i++) 
     stream<<(*ii)[i]<<" ";
     stream<<'\n';
     }
     return stream;
     }
     */
};

#endif

