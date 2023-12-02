/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! \file AddQQVspJT.cxx
*/
#ifndef __ADDQQVSPJT__CXX__
#define __ADDQQVSPJT__CXX__

/* ------- This file only exists because putting both Particle - Particle and Particle - Hole
 Formalisms in the same file (spherical2sp.cxx) seems like a not very optimal idea.
 Furthermore, There are some changes that may need to be made so not editing original files
 might prove redeeming in the future.
 
 This file dublicates AddQQVsp
 
 */

//#include <tav/tensor.cxx>
#include <tav/stensor.cxx>
#include <tav/ThreeJSymbolM.cxx>
using tav::ThreeJSymbol;
#include <sm/QN.cxx>
#include <sm/SingleParticleStates.cxx>
//#include <sm/MBSQ.cxx>
#include "cioperators.cxx"
#include <tav/pconfiguration.cxx>
//using SM::SPSqrM;
namespace SM{
    
    int AddQQVspJT(
                   tav::STensor<spsint,double> &ESP,  ///< Output map of e[12] elements
                   tav::STensor<spsint,double> &VSP,  ///< Output map of elements V[1234]  (12) (34)-ordered pairs
                   SingleParticleStates & sps, ///< s.p. states matrix
                   const int *ijk, ///< list 1 2 3 4 L T
                   const double V,  ///< interaction strength
                   bool usepn = false,//set weather p-n will be specifically selected
                   int pnstate = 0//controls Tz projection allowed
    )
    {
        double htmp,cgj1,cgj2,cgt1,cgt2,cgJT;
        int bl;
        bool phase;
        spsint a[4];
        spsint b[4];
//        cgJ = mav::ClebschGordan(
        //order of labels 0123 but operators
        //a0+ a1 a3+ a2
        //here we implement a brute force scan
        for (a[0]=0;a[0]<sps.size();a[0]++) {
            if (sps[a[0]].level!=ijk[0]) continue;
            for (a[1]=0;a[1]<sps.size();a[1]++) {
                if (sps[a[1]].level!=ijk[1]) continue;
                for (a[2]=0;a[2]<sps.size();a[2]++) {
                    if (sps[a[2]].level!=ijk[2]) continue;
                    for (a[3]=0;a[3]<sps.size();a[3]++) {
                        if (sps[a[3]].level!=ijk[3]) continue;
                        //check z projections
                        if (sps[a[0]].Jz+sps[a[3]].Jz!=sps[a[1]].Jz+sps[a[2]].Jz) continue;
                        if (sps[a[0]].Tz+sps[a[3]].Tz!=sps[a[1]].Tz+sps[a[2]].Tz) continue;
                        
                        if (usepn && int(sps[a[0]].Tz+sps[a[3]].Tz) != 2*pnstate) continue; //this is extra condition on pn state to select only certain typse of pairs
                        //spin
                        cgj1 = mav::ClebschGordan(sps[a[0]].J/2., sps[a[0]].Jz/2.,
                                                  sps[a[1]].J/2., -sps[a[1]].Jz/2.,
                                                  ijk[4], (sps[a[0]].Jz - sps[a[1]].Jz)/2.);
                        cgj2 = mav::ClebschGordan(sps[a[2]].J/2., sps[a[2]].Jz/2.,
                                                  sps[a[3]].J/2., -sps[a[3]].Jz/2.,
                                                  ijk[4], (sps[a[2]].Jz - sps[a[3]].Jz)/2.);
                        if (std::abs(cgj1) < 1.e-10 || std::abs(cgj2) < 1.e-10) continue;
                        cgt1 = mav::ClebschGordan(sps[a[0]].T/2., sps[a[0]].Tz/2.,
                                                  sps[a[1]].T/2., -sps[a[1]].Tz/2.,
                                                  ijk[5], (sps[a[0]].Tz - sps[a[1]].Tz)/2.);
                        cgt2 = mav::ClebschGordan(sps[a[2]].T/2., sps[a[2]].Tz/2.,
                                                  sps[a[3]].T/2., -sps[a[3]].Tz/2.,
                                                  ijk[5], (sps[a[2]].Tz - sps[a[3]].Tz)/2.);
                        if (std::abs(cgt1) < 1.e-10 || std::abs(cgt2) < 1.e-10) continue;
                        
                        
                        //here phase is used just for the first 2 clebsches to not call the function.
                        phase = (((2 * (ijk[4] + ijk[5]) - (sps[a[0]].Jz - sps[a[1]].Jz + sps[a[0]].Tz - sps[a[1]].Tz))>>1)&1);//(-)^{J - μ + Τ - τ}
                        cgJT = 1./(std::sqrt((2 * ijk[4] + 1.) * (2 * ijk[5] + 1.)));
                        cgJT = ( phase ? -cgJT : cgJT );
//                        cgJT = mav::ClebschGordan(ijk[4],(sps[a[0]].Jz - sps[a[1]].Jz)/2.,
//                                                  ijk[4],(sps[a[2]].Jz - sps[a[3]].Jz)/2.,
//                                                  0,0);
//                        cgJT*= mav::ClebschGordan(ijk[5],(sps[a[0]].Tz - sps[a[1]].Tz)/2.,
//                                                  ijk[5],(sps[a[2]].Tz - sps[a[3]].Tz)/2.,
//                                                  0,0);
                        htmp = cgJT * cgj1 * cgj2 * cgt1 * cgt2;
                        if (std::abs(htmp) < 1.e-10) continue;

                        phase  = ((sps[a[1]].J + sps[a[1]].Jz + sps[a[1]].T + sps[a[1]].Tz)>>1)&1;
                        phase ^= ((sps[a[3]].J + sps[a[3]].Jz + sps[a[3]].T + sps[a[3]].Tz)>>1)&1;
                        //if
                        if (a[1]==a[2])
                            ESP(a[0],a[3]) += (phase ? -htmp * V : htmp * V);
                        
                        //Not allowed by Fermi Principle
                        if (a[0] == a[2]) continue;
                        if (a[1] == a[3]) continue;
                        b[0] = a[0]; b[1] = a[2]; b[2] = a[1]; b[3] = a[3];
                        if(b[0]>b[1]) {std::swap(b[0],b[1]); phase=!phase;}
                        if(b[2]>b[3]) {std::swap(b[2],b[3]); phase=!phase;}
                        if(b[0] > b[2] && b[1] > b[3])
                        {
                            std::swap(b[0],b[2]);
                            std::swap(b[1],b[3]);
                        }
                        //There is an extra (-) sign forom commutation.
                        VSP[b]+=(phase ? htmp * V : -htmp * V);

                    }//end a[3]
                }//end a[2]
            }//end a[1]
        }//end a[0]
        StripZeros(VSP);
        return 0;
    }//end function
    
}//end namespace SM


#endif //__ADDQQVSPJT__CXX__
