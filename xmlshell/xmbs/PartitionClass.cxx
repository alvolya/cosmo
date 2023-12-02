/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


#ifndef __PARTITION__CLASS__CXX__
#define __PARTITION__CLASS__CXX__

#include <list>
#include <vector>
#include <iostream>
using namespace std;

#define FatalError_PartitionClass(message)         { std::cerr <<"Error: "<< message <<std::endl<< __FILE__ << ":" << __LINE__ << std::endl; throw (1); }

/** This code checks max-min occupation rejections
 true=reject, false-continue
 similar to Knapsack problem?
 4/19/20 remove dependence on Min
 */
bool TestRestriction(
                     uschar *a, ///<test array
                     std::vector< std::vector<bool> > &rej, ///< list of rejections
                     std::vector<int> &min, ///<list of min
                     std::vector<int> &max ///list of max
)
{
    
    int sum=0;
    for (int rr=0;rr<min.size();rr++)
    {
        sum=0;
        for (int i=0;i<rej[rr].size();i++) if (rej[rr][i]) sum+=a[i];
        if ((sum<min[rr])||(sum>max[rr])) {return true;}
    }
    return false;
}




/** This is a class that generates partitions
 */
class Partition: public vector<uschar *>
{
private:
    int *partition_;
    void Print_(
                int N, ///<particle number
                int lev ///<number of levels
    )
    {
        //cout<<"Enter with: "<<N<<" "<<lev<<endl;
        if (lev==1) // one level simple
        {
            if ((N<=wmax[0])&&(N>=wmin[0]))
            {
                partition_[0]=N;
                //once we get here we have a good partition
                for (int i=0;i<parts;i++) cout<<partition_[i]<<" ";
                cout<<endl;
            }
        }
        else {
            for (int nn=wmin[lev-1];((nn<=wmax[lev-1])&&(nn<=N));nn++)
            {
                //cout<<nn<<endl;
                partition_[lev-1]=nn;
                Print_(N-nn,lev-1);
            }
        }
    };
    /**  Make partition priority right ascending,
     */
    void Make_(
               int N, ///<particle number
               int lev ///<number of levels
    )
    {
        // cerr<<" make "<<N<<endl;
//        ShowStatus();
        if (lev==1) // one level simple
        {
            if ((N<=wmax[0])&&(N>=wmin[0])) {
                partition_[0]=N;
                //once we get here we have a good partition
                // for (int i=0;i<parts;i++) cout<<partition_[i]<<" ";
                // cout<<endl;
                // push_back(partition_);
                //partition_=new int [parts];
                uschar *tmp =new uschar [parts];
                for (int q=0;q<parts;q++) tmp[q]=partition_[q];
                push_back(tmp);
                
            }
        }
        else {
            for (int nn=wmin[lev-1];((nn<=wmax[lev-1])&&(nn<=N));nn++){
                partition_[lev-1]=nn;
                Make_(N-nn,lev-1);
            }
        }
    };

    // This is unnecessary because makes Qmin=0 and Qmax=Nmax+Ncore
//    //dependence on Nmax->Nmax+Ncore
//    void Make_(
//                int N,///<particle number
//                int lev,///<number of levels
//                int Nmax,///<maximum number of allowed quanta;
////                int NCore,///<quanta for lowest configuration
//                vector<int> qq///<vector of size lev with the quanta for each level
//    )
//    {
//        if (lev==1) // one level simple
//        {
//            if ((N<=wmax[0])&&(N>=wmin[0]))
//            {
//                partition_[0]=N;
//                //once we get here we have a good partition
//                // for (int i=0;i<parts;i++) cout<<partition_[i]<<" ";
//                // cout<<endl;
//                // push_back(partition_);
//                //partition_=new int [parts];
//                uschar *tmp =new uschar [parts];
//                int qqn=0;
//                for (int q=0;q<parts;q++)
//                {
//                    tmp[q]=partition_[q];
////                    cout<<partition_[q]<<" ";
//                    qqn+=partition_[q]*qq[q];
//                }
////                cout<<"("<<qqn<<")"<<endl;
//                if (qqn<=Nmax)
//                    push_back(tmp);
//
//            }
//        }
//        else {
//            for (int nn=wmin[lev-1];((nn<=wmax[lev-1])&&(nn<=N));nn++){
//                partition_[lev-1]=nn;
//                Make_(N-nn,lev-1,Nmax,qq);
//            }
//        }
//    }
    
    /*This makes partition with limitations on Abelian quanta
     */
        void MakeQ_(
                    int N,///<particle number
                    int lev,///<number of levels
                    int Qmin,///<minimum number of quanta
                    int Qmax,///<maximum number of quanta
                    const vector<int> &qq///<vector of size lev with the quanta for each level
        )
        {
            //cerr<<"Call: "<<N<<" "<<lev<<" "<<Qmin<<" "<<Qmax<<endl;
            if (lev==1) // one level simple
            {
                int qqn=N*qq[0];
                if ((N<=wmax[0])&&(N>=wmin[0])&&(qqn<=Qmax)&&(qqn>=Qmin))
                {
                    partition_[0]=N;
                    //once we get here we have a good partition
                    // for (int i=0;i<parts;i++) cout<<partition_[i]<<" ";
                    // cout<<endl;
                    // push_back(partition_);
                    //partition_=new int [parts];
                    uschar *tmp =new uschar [parts];
                    for (int q=0;q<parts;q++)
                    {
                        tmp[q]=partition_[q];
                    }
                        push_back(tmp);
                        
                }
            }
            else {
                for (int nn=wmin[lev-1];((nn<=wmax[lev-1])&&(nn<=N));nn++){
                    //put nn particles on last level
                    int qremove=nn*qq[lev-1];
                    partition_[lev-1]=nn;
                    MakeQ_(N-nn,lev-1,Qmin-qremove,Qmax-qremove,qq);
                }
            }
        }
    
    
    /**  Make partition priority right ascending,
     \todo there is some kind of bug with make2, does not work for spin projections
     */
    void Make2_(
                int N, ///<particle number
                int lev, ///<number of levels
                int *px,
                int *pmx,
                int *pmn ///<set of pointers
    )
    {
        if (lev==1) // one level simple
        {
            if ((N<=(*pmx))&&(N>=(*pmn))) {
                (*px)=N;
                //once we get here we have a good partition
                uschar *tmp =new uschar [parts];
                for (int q=0;q<parts;q++) tmp[q]=partition_[q];
                push_back(tmp);
                //for (int i=0;i<parts;i++) cout<<partition_[i]<<" ";
                //cout<<endl;
            }
        }
        else {
            //for (int nn=pmn[0];((nn<=pmx[0])&&(nn<=N));nn++){
            //for (int nn=tav::Min(*wmax,N); nn>=(*wmin);nn--){
            for (int nn=((*wmax<N)?(*wmax):N) ; nn>=(*wmin);nn--){
                (*px)=nn;
                Make2_(N-nn,lev-1,px+1,pmx+1,pmn+1);
            }
        }
    };
    
public:
    int sum;
    int *wmax;
    int *wmin;
    int parts;
    int Initialize(int s,int p, const int *wx, const int *wm) {
        sum=s; parts=p;
        wmax=new int [parts];
        wmin=new int [parts];
        partition_= new int [parts];
        for (int i=0;i<parts;i++) {wmax[i]=wx[i]; wmin[i]=wm[i];}
        return 0;
    }
    
    Partition() {
        sum=0; parts=0;
    }
    
    Partition(int s,int p) {
        sum=s; parts=p;
        wmax=new int [parts];
        wmin=new int [parts];
        partition_= new int [parts];
        for (int i=0;i<parts;i++) {wmax[i]=sum;wmin[i]=0;}
    }
    
    Partition(int s,int p, int *wx) {
        sum=s; parts=p;
        wmax=new int [parts];
        wmin=new int [parts];
        partition_= new int [parts];
        for (int i=0;i<parts;i++) {wmax[i]=wx[i]; wmin[i]=0;}
    }
    
    Partition(int s,int p, const int *wx, const int *wm) {
        sum=s; parts=p;
        wmax=new int [parts];
        wmin=new int [parts];
        partition_= new int [parts];
        for (int i=0;i<parts;i++) {wmax[i]=wx[i]; wmin[i]=wm[i];}
    }
    
    ~Partition(void){
        delete [] wmax;
        delete [] wmin;
        delete [] partition_;
        iterator ii;
        for (ii=begin(); ii!=end();ii++) delete [] *ii;
    }
    /** Insert into my partition other partitions P and P
     */
    Partition(Partition &P, Partition &Q)
    {
        sum=P.sum+Q.sum;
        parts=P.parts+Q.parts;
        wmax=new int [parts];
        wmin=new int [parts];
        partition_= new int [parts];
        // cout<<parts<<" "<<P.size()<<" "<<Q.size()<<endl;
        //Pause();
        int i;
        for (i=0;i<P.parts;i++)
        {
            wmax[i]=P.wmax[i];
            wmin[i]=P.wmax[i];
        }
        for (;i<parts;i++)
        {
            wmax[i]=Q.wmax[i];
            wmin[i]=Q.wmin[i];
        }
        for (Partition::iterator ip=P.begin(); ip!=P.end();ip++)
            for (Partition::iterator iq=Q.begin(); iq!=Q.end();iq++)
            {
                uschar *tmp;
                try
                {
                    tmp =new uschar [parts];
                }
                catch (std::exception& e)
                {
                    cerr<<"Exception: Partition(Partition P, Partitition Q);"<<endl;
                    FatalError_PartitionClass(e.what());
                }
                
                for (i=0;i<P.parts;i++) tmp[i]= (*ip)[i];
                for (int q=0;q<Q.parts;q++) {tmp[i]= (*iq)[q]; i++;}
                push_back(tmp);
            }
    }
    
    
    
    void Print(void) {Print_(sum,parts);}
    int Make(void) {Make_(sum,parts); return size();}
    //int Make(int Nmax,vector<int> qq) {Make_(sum,parts,Nmax,qq); return size();}
    int Make(int Nmax,vector<int> qq) {MakeQ_(sum,parts,0,Nmax,qq); return size();}
    int MakeQ(int Qmin,int Qmax, const vector<int> &qq) {MakeQ_(sum,parts,Qmin,Qmax,qq); return size();}
    int Make2(void) {Make2_(sum,parts,partition_,wmax,wmin); return size();}
    
    friend std::ostream& operator << (std::ostream &stream,Partition &P)
    {
        for (Partition::iterator ii=P.begin(); ii!=P.end();ii++)
        {
            for (int i=0;i<P.parts;i++)
                stream<<int((*ii)[i])<<" ";
            stream<<'\n';
        }
        return stream;
    }
};
/** Place parity restriction of a given partition
 */
int ParityPurge(Partition &P,///<Partition to purge
                int no,///<number of levels of negative parity
                const int *orbit, ///<list of levels of negative parity
                bool parity ///<parity that we want (true=negative)
){
    bool myparity=false;//positive
    for (Partition::iterator ii=P.begin(); ii!=P.end();ii++)
    { myparity=false; //positive
        for (int i=0;i<no;i++) myparity^=((*ii)[orbit[i]]&1);
        if (myparity^parity) {delete [] (*ii); ii=P.erase(ii); ii--;};
    }
    return P.size();
}
/** Place parity restriction of a given partition, list version
 Parity \f$ (p)^n \f$ is logical AND operation p&&(n&1)
 */
int ParityPurge(Partition &P,///<Partition to purge
                bool *orbit, ///<list of levels of negative parity
                bool parity ///<parity that we want (true=negative)
){
    bool myparity=false;//positive
    for (Partition::iterator ii=P.begin(); ii!=P.end();ii++)
    { myparity=parity; //positive
        for (int i=0;i<P.parts;i++) myparity^=orbit[i]&&(((*ii)[i])&1);
        if (myparity) {delete [] (*ii); ii=P.erase(ii); ii--;};
    }
    return P.size();
}
//This is unnecessary because use min=0 max=Nmax+Ncore
//int QuantaPurge(Partition &P,///Partition to purge
//                vector<int> &quanta,///list of quanta per level
//                int Nmax,///Allowed number of quanta
//                int NCore///The number of core quanta
//){
//    int qq=0;
//    for (Partition::iterator it=P.begin();it!=P.end();it++)
//    {
//        qq=0;
//        for (int i=0;i<P.parts;i++)
//        {
//            qq+=quanta[i]*int((*it)[i]);
//        }
////        cout<<qq<<" "<<NCore<<" "<<Nmax<<endl;
//        if (qq-NCore>Nmax)
//        {
//            delete [] (*it);
//            it= P.erase(it);
//            it--;
//        }
//
//    }
//    return P.size();
//}
int QuantaMinMaxPurge(Partition &P,///Partition to purge
                vector<int> &quanta,///list of quanta per level
                int min,///Allowed number of quanta
                int max///The number of core quanta
){
    int qq=0;
    for (Partition::iterator it=P.begin();it!=P.end();it++)
    {
        qq=0;
        for (int i=0;i<P.parts;i++)
        {
            qq+=quanta[i]*int((*it)[i]);
        }
        //        cout<<qq<<" "<<NCore<<" "<<Nmax<<endl;
        if (qq>max || qq<min)
        {
            delete [] (*it);
            it= P.erase(it);
            it--;
        }
        
    }
    return P.size();
}
/** Place restriciton on a partition
 */
int RestrictionPurge(Partition &P,///<Partition to purge
                     int no,///<number of levels
                     const int *orbit, ///<list of levels
                     int min, ///<Minimum of the sum
                     int max ///< Maximum of the sum
){
    for (Partition::iterator ii=P.begin(); ii!=P.end();ii++)
    { int sum=0;
        for (int i=0;i<no;i++) sum+=((*ii)[orbit[i]]);
        if ((sum<min)||(sum>max)) {delete [] (*ii); ii=P.erase(ii); ii--;};
    }
    return P.size();
}

/** Place restriciton on a partition with boolean list
 */
int RestrictionPurge(Partition &P,///<Partition to purge
                     bool *orbit, ///<list of levels
                     int min, ///<Minimum of the sum
                     int max ///< Maximum of the sum
){
    for (Partition::iterator ii=P.begin(); ii!=P.end();ii++)
    { int sum=0;
        for (int i=0;i<P.parts;i++) if (orbit[i]) sum+=((*ii)[i]);
        if ((sum<min)||(sum>max)) {delete [] (*ii); ii=P.erase(ii); ii--;};
    }
    return P.size();
}


int MakePartition(Partition &Z,  Partition &P, Partition &Q)
{
    
    int sum=P.sum+Q.sum;
    int parts=P.parts+Q.parts;
//    cout<<Q.sum<<" "<<Q.parts<<endl;
    int* wmax=new int [parts];
    int* wmin=new int [parts];
    // cout<<parts<<" "<<P.size()<<" "<<Q.size()<<endl;

    //Pause();
    int i;
    for (i=0;i<P.parts;i++) {wmax[i]=P.wmax[i]; wmin[i]=P.wmin[i];}
    for (;i<parts;i++) {wmax[i]=Q.wmax[i]; wmin[i]=Q.wmin[i];}

    Z.Initialize(sum,parts,wmin,wmax);
    
    for (Partition::iterator ip=P.begin(); ip!=P.end();ip++)
    {
        for (Partition::iterator iq=Q.begin(); iq!=Q.end();iq++)
        {
            uschar *tmp;
            try {tmp =new uschar [parts];}
            catch (std::exception& e)
            {
                cerr<<"Exception: MakePartition(Partition P, Partitition Q);"<<endl;
                FatalError_PartitionClass(e.what());
            }
            
            for (i=0;i<P.parts;i++) tmp[i]= (*ip)[i];
            for (int q=0;q<Q.parts;q++) {tmp[i]= (*iq)[q]; i++;}
            Z.push_back(tmp);
        }
    }

    delete [] wmax;
    delete [] wmin;
    return 0;
}
/** Combine two partitions and use rejections
 */
int MakePartition(Partition &Z,///<new partition
                  Partition &P, ///<1st partition
                  Partition &Q, ///<2nd partition
                  std::vector< std::vector<bool> > &rej, ///< list of rejections
                  std::vector<int> &min, ///<list of min
                  std::vector<int> &max ///<list of max
)
{
    
    int sum=P.sum+Q.sum;
    int parts=P.parts+Q.parts;
    int* wmax=new int [parts];
    int* wmin=new int [parts];
    // cout<<parts<<" "<<P.size()<<" "<<Q.size()<<endl;
    //Pause();
    int i;
    for (i=0;i<P.parts;i++) {wmax[i]=P.wmax[i]; wmin[i]=P.wmin[i];}
    for (;i<parts;i++) {wmax[i]=Q.wmax[i]; wmin[i]=Q.wmin[i];}
    
    Z.Initialize(sum,parts,wmin,wmax);
    
    for (Partition::iterator ip=P.begin(); ip!=P.end();ip++)
    {
        for (Partition::iterator iq=Q.begin(); iq!=Q.end();iq++)
        {
            //ShowStatus();
            uschar *tmp;
            try {tmp =new uschar [parts];}
            catch (std::exception& e)
            {
                cerr<<"Exception: MakePartition(Partition P, Partitition Q, rejctions);"<<endl;
                FatalError_PartitionClass(e.what());
            }
            
            for (i=0;i<P.parts;i++) tmp[i]= (*ip)[i];
            for (int q=0;q<Q.parts;q++) {tmp[i]= (*iq)[q]; i++;}
            if(TestRestriction( tmp, rej, min, max )) {delete [] tmp;} 
            else Z.push_back(tmp);
        }
    }
    delete [] wmax;
    delete [] wmin;
    return 0;
}

/** Combine two partitions and use rejections
 */
int MakePartitionOMP(Partition &Z,///<new partition
                  Partition &P, ///<1st partition
                  Partition &Q, ///<2nd partition
                  std::vector< std::vector<bool> > &rej, ///< list of rejections
                  std::vector<int> &min, ///<list of min
                  std::vector<int> &max ///<list of max
)
{
//    cout<<"In here"<<endl;
    int sum=P.sum+Q.sum;
    int parts=P.parts+Q.parts;
    int* wmax=new int [parts];
    int* wmin=new int [parts];
    // cout<<parts<<" "<<P.size()<<" "<<Q.size()<<endl;
    //Pause();
    int i;
    for (i=0;i<P.parts;i++) {wmax[i]=P.wmax[i]; wmin[i]=P.wmin[i];}
    i=0;
    for (;i<parts;i++) {wmax[i]=Q.wmax[i]; wmin[i]=Q.wmin[i];}
    
    Z.Initialize(sum,parts,wmin,wmax);
    //Make vectors of iterators.
    vector<Partition::iterator> ip_vec(P.size()),iq_vec(Q.size());
    i=0;
    for (Partition::iterator itq=P.begin();itq!=P.end();itq++,i++) ip_vec[i]=itq;
    i=0;
    for (Partition::iterator itq=Q.begin();itq!=Q.end();itq++,i++) iq_vec[i]=itq;
    //first pick which partition is bigger to parallelize over.
    
    std::vector<std::list<std::vector<uschar> > > Z_list(P.size()>Q.size() ? ip_vec.size() : iq_vec.size());
    cout<<Q.size()<<" - "<<P.size()<<endl;

    if (P.size()>Q.size())//parallelize P.
    {
#pragma omp parallel for
        for (int j=0;j<ip_vec.size();j++)
        {
            int ii=0;
            vector<uschar> tmp(parts);
            Partition::iterator ip=ip_vec[j];
            for (int k=0;k<iq_vec.size();k++)
            {
                Partition::iterator iq=iq_vec[k];
                for (ii=0;ii<P.parts;ii++) tmp[ii]= (*ip)[ii];
                for (int q=0;q<Q.parts;q++,ii++) {tmp[ii]= (*iq)[q];}
                bool truefalse=TestRestriction( &(tmp[0]), rej, min, max );
                if (!truefalse) Z_list[j].push_back(tmp);

            }
        }
        
    }
    else//Parallelize Q.
    {
#pragma omp parallel for
        for (int j=0;j<iq_vec.size();j++)
        {
            int ii=0;
            vector<uschar> tmp(parts);
            Partition::iterator iq=iq_vec[j];
            for (int k=0;k<ip_vec.size();k++)
            {
                Partition::iterator ip=ip_vec[k];
                for (ii=0;ii<P.parts;ii++) tmp[ii]= (*ip)[ii];
                for (int q=0;q<Q.parts;q++,ii++) {tmp[ii]= (*iq)[q];}
                bool truefalse=TestRestriction( &(tmp[0]), rej, min, max );
                if (!truefalse) Z_list[j].push_back(tmp);
            }
        }
    }

    //And now to safely copy back everything.
    for (int q=0;q<Z_list.size();q++)
    {
        std::list<vector<uschar> >::iterator itl=Z_list[q].begin();
        for (;itl!=Z_list[q].end();itl++)
        {
            uschar *tmp;
            try {tmp =new uschar [parts];}
            catch (std::exception& e)
            {
                cerr<<"Exception: MakePartition(Partition P, Partitition Q, rejctions);"<<endl;
                FatalError_PartitionClass(e.what());
            }
            for (int qq=0;qq<parts;qq++)
                tmp[qq]=(*itl)[qq];
            Z.push_back(tmp);
        }
    }
//    cout<<Z_list.size()<<" partitions\n";
    delete [] wmax;
    delete [] wmin;
    return 0;
}
#endif



