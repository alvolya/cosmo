/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*
 06/30/2008 Bug fix in beta decay
 this is modification of, still needs work ... EMSPMul.cxx
 is program to find magnetic and electric single-particle multipoles
 07/21/2006 I have decided to use BM definition of operator M see page 384
 Eq. 3C-36
 output is MUM[mf][mi]=<f|M|i>, note there is no i^\lambda
 magnetic operator is in units of nuclear magneton e\hbar/(2Mc)
 For magnetic moment see Eq. 3-40 page 337 Bohr Mottelson
 tested against data  Phys. Rev. C62, 044312 (2000).
 BGT have passed a single-particle test BM eq. 3-47 page 346
 06/10/2016 Bug fix: initial and final states were flipped in ROperator_YS
 Verification tested in sd, see B. H. Wildenthal, M. S. Curtin, and B. A. Brown, Phys. Rev. C. Phys. Rev. C 28, 1343 (1983)
 
 */
/**
 \author Alexander Volya  <http://www.volya.net>
 \file EMspherical2sp.cxx
 \brief Electromagnetic and weak operators
 \example SHLBEM.cpp
 */
#ifndef __EMSPHERICAL2SP__CXX__
#define __EMSPHERICAL2SP__CXX__
#ifndef PI
#define PI 3.141592653589793238462643
#endif
#include <constants.cxx>
#include <tav/tensor.cxx>
#include <sm/cioperators.cxx> //use for stensor class
#include <mav/SixJSymbol.cxx>
using mav::SixJSymbol;
#include <tav/ThreeJSymbolM.cxx> //threeJ and clebsch are required
//#include <mav/ThreeJSymbol.cxx> //threeJ and clebsch are required
using tav::ClebschGordan;
//#include <sm/System.cxx> //system is needed
#include <sm/SingleParticleStates.cxx> //we need single particle states
#include <sm/HO.cxx> //harmonic oscillator radii
#ifndef PARITY
#define PARITY(a) (1-(((a)&1)<<1))
#endif
namespace SM{
    /** Given reduced matrix element compute full matrix element
     \f[ \langle j'm'|Q_{\lambda\mu}^{\dagger}|jm\rangle= \frac{(-)^{2\lambda}}{\sqrt{2j'+1}}C_{jm\lambda\mu}^{j'm'}\, \langle j'||Q_{\lambda}^{\dagger}||j\rangle \f]
     Non-zero M is only computed
     
     This form is consistent with Edmonds (page 75). However, full matrix element in Bohr-Mottelson (eq 1A-60, page 82) has no extra phase
     \f$ (-1)^{2\lambda} \f$ which has no effect on integers-spin operators.
     Below is a set of related definitions.
     \f[
     \langle \frac{\langle j'm'|Q_{\lambda\mu}^{\dagger}|jm\rangle}{\langle j'||Q_{\lambda}^{\dagger}||j\rangle}=(-)^{j'-m'}\left(\begin{array}{ccc}
     j' & \lambda & j\\
     -m' & \mu & m\end{array}\right)=\frac{(-)^{j-m}}{\sqrt{2\lambda+1}}C_{j'm'j-m}^{\lambda\mu}=\frac{(-)^{\lambda-j+j'}}{\sqrt{2j'+1}}C_{\lambda\mu jm}^{j'm'}=\frac{(-)^{2\lambda}}{\sqrt{2j'+1}}C_{jm\lambda\mu}^{j'm'}
     \f]
     */
    double ReducedToFull(
                         QN &JF, ///< Final sef of quantum numbers
                         QN &JI,  ///< Initial set of quantum numbers
                         int Lambda, ///< Multipole Lambda-integer, not 2*spin format
                         double reducedme=1.0 ///<Reduced matrix element
    )
    {
        return reducedme*tav::ClebschGordan(JI.J, JI.Jz,
                                            2*Lambda, (JF.Jz-JI.Jz),
                                            JF.J, JF.Jz
                                            )/sqrt(double(JF.J+1));
        
    }
    
    
    
    /** \brief Magnetic multipole
     \f[
     \langle 2||{\cal M}(M\lambda)||1 \rangle = \mu_N \, i^{l_1-l_2}\,(-1)^{j_1-j_2+\lambda-1}
     \sqrt{\frac{(2\lambda+1)(2j_1+1)}{4\pi}} \langle 2|r^{\lambda-1}|1\rangle  (A+B)
     \f]
     \f[
     A=\left ( g_s -\frac{2g_l}{\lambda+1}\right )\frac{\lambda}{2}\, \langle j_1 \frac{1}{2}, \lambda 0 |j_2 \frac{1}{2} \rangle \left [ 1+\frac{1}{\lambda} (-1)^{l_1+1/2-j_1} \left \{(j_1+1/2)+(-1)^{j_1+j_2-\lambda}(j_2+1/2)\right \}\right ]
     \f]
     \f[
     B=(-1)^{j_1+j_2+\lambda} \frac{2g_l}{\lambda+1} \sqrt{\lambda(4\lambda^2-1)j_1(j_1+1)(2j_1+1)}\,\, \langle j_1 \frac{1}{2}, \lambda-1 0 |j_2 \frac{1}{2} \rangle \,\,
     \left\{ \begin{array}{ccc}
     j_{1} & 1 & j_{1}\\
     \lambda-1 & j_{2} & \lambda
     \end{array}\right\}
     \f]
     
     Matrix element is in units of nuclear magneton  \f$ \mu_N=\frac {e \hbar}{2 M_p c}\f$
     magnetic moments are expected in the same units; the  g factor is twice
     that, because spin is 1/2 \mu=g s. see Particle Data Group for data. \f$ g_l=0 \f$ for neutron and \f$ g_l=1 \f$ for proton
     See also Bohr-Mottelson vol I, Eq. 3C-36, page 388.
     
     The magnetic moment (vector) operator is \f[\mu_z=\sqrt{\frac{4\pi}{3}} {\cal M}(M1)_{\lambda \mu=0} \f]
     
     
     */
    template <class Radial>
    int MagneticSPMultipole(
                            MBOperator &MUM, ///<output operator
                            SingleParticleStates &sps, ///<single particle states
                            int K, ///< Momentum of the operator
                            int kappa, ///< Magnetic projection of the operator
                            Radial &Rfi, ///< Matrix of radial overlaps \f$ r^{K-1} \f$
                            double mun=phy::muneutron, ///< Magnetic moment of the neutron
                            double mup=phy::muproton, ///< Magnetic moment of the proton
                            double GL_neutron=0.0,  /// Gl for neutron, use defult unles you need something else
                            double GL_proton=1.0 /// Gl value for proton
    ){
        spsint pp[2]; //we first do pairs, must use constants
        /*sps -single particle states
         Rfi -radial overlap integral <f|r^{K-1}|i> convention R(r)>0 r-> \infty
         mun -effective neutron magnetic moment in terms of nuclear magneton
         mup -proton matnetic moment in nuclear magnetons
         nuclear magneton = e\hbar/(2 M_p c)
         magnetic moments are also in the same units, however g factor is twice
         that, because spin is 1/2 \mu=g s. see Particle Data Group for data.
         Bohr-Mottelson for definition and units
         K-multipolatity
         kappa-projection of multipolatiry
         output Matirix in SP basis
         */
        double GS,GL;
        double htmp, htmp1,htmp2;
        for (int mf=0;mf<sps.size();mf++)
            for (int mi=0;mi<sps.size();mi++) {
                //check parity and other simple stuff
                if  (sps[mf].Tz!=sps[mi].Tz) continue;
                if (((sps[mi].L-sps[mf].L+K-1)%2!=0)) continue;
                if ((sps[mi].Jz+2*kappa-sps[mf].Jz)!=0) continue;
                //reduced matrix element, page 387, volume I
                // if ((sps[mf].Tz)==1) {GL=0.0; GS=2.0*phy::muneutron;}
                // else {GS=2.0*phy::muproton; GL=1.0;}
                if ((sps[mf].Tz)==1) {GL=GL_neutron; GS=2.0*mun;}
                else {GS=2.0*mup; GL=GL_proton;}
                htmp1=PARITY((((sps[mi].J+sps[mf].J)>>1)-K))*(sps[mf].J+1.0)/2.0+
                (sps[mi].J+1.0)/2.0;
                htmp1*=PARITY(sps[mi].L-(sps[mi].J-1)/2)/double(K);
                htmp1+=1.0;
                htmp1*=tav::ClebschGordan(sps[mi].J, 1,
                                          (2*K), 0,
                                          sps[mf].J, 1 );
                htmp1*=double(K)/2.0;
                htmp1*=(GS-2.0*GL/(double(K)+1));
                //---------------------end of first term---------
                htmp2=SixJSymbol(sps[mi].J/2., 1.0, sps[mi].J/2.,
                                 double(K)-1.0, sps[mf].J/2.0, double(K));
                htmp2*= tav::ClebschGordan(sps[mi].J, 1,
                                           2*K-2, 0,
                                           sps[mf].J, 1 );
                htmp2*=sqrt(K*(4.0*K*K-1)*sps[mi].J*(sps[mi].J+2.0)*(sps[mi].J+1))*GL;
                htmp2*=PARITY((((sps[mi].J+sps[mf].J)>>1)+K))/(double(K)+1.0);
                //---------------------end second term------
                htmp=htmp1+htmp2;
                if (((((sps[mi].J-sps[mf].J)>>1)+K-1)%2)!=0) htmp=-htmp;
                //this is phase factor from i^{li-lf+K-1} removed now.
//                if (((sps[mi].L-sps[mf].L)%4)!=0) htmp=-htmp;
                htmp*=sqrt((2.0*K+1.)*(sps[mi].J+1.)/(4.0*PI));
                htmp*=Rfi[sps[mf].level][sps[mi].level]; //now I multiply by Rfi
                //if (std::abs(htmp)<1E-10)  continue; //zero due to radial happens often
                if (htmp==0)  continue;
                //cout<<htmp<<endl;
                
                
                //reduced matrix element is done
                //full matrix element, eq 1A-60, page 82
                /*
                 htmp*= ClebschGordan(sps[mi].J/2., sps[mi].Jz/2.,
                 double(K), double(kappa),
                 sps[mf].J/2., sps[mf].Jz/2.
                 );
                 
                 htmp/=sqrt(double(sps[mf].J+1));
                 */
                //   if ((mi==mf)&&(sps[mi].J==sps[mi].Jz)) cerr<<sps[mi].Tz<<" "<<sps[mi].J/2.0<<" "<<htmp*sqrt(4.0*PI/3.0)<<endl; //for M1 these are sp moments
                pp[0]=mi; pp[1]=mf;
                MUM[pp]+=ReducedToFull(sps[mf],sps[mi],K,htmp);
            }//end loop
        
        return 0;
    }
    
    ///\brief Magnetic multipole with Harmonic Oscillator wf.
    /*
     (1/5/2016)--> Added a phase convention bool variable which is true by default but when set to true assumes wavefunctions that are positive at infinity (extra factor of \f$ (-1)^{n_1 + n_2} \f$
     */
    
    int MagneticSPMultipole(
                            MBOperator &MUM, ///<output operator
                            SingleParticleStates &sps, ///<single particle states
                            int K, ///< Momentum of the operator
                            int kappa, ///< Magnetic projection of the operator
                            double mun=phy::muneutron, ///< Magnetic moment of the neutron
                            double mup=phy::muproton, ///< Magnetic moment of the proton
                            double GL_neutron=0.0,  /// Gl for neutron, use defult unles you need something else
                            double GL_proton=1.0, /// Gl value for proton
                            bool PhaseConvention=true) {
        // Determine the number of levels
        int levels=0;
        for (int mf=0;mf<sps.size();mf++) if (levels<sps[mf].level) levels=sps[mf].level;
        levels++;
        // create matrix
        tav::Tensor<2,double> Rfi(levels,levels);
        for (int mf=0;mf<sps.size();mf++)
            for (int mi=0;mi<sps.size();mi++)
                Rfi[sps[mf].level][sps[mi].level]=
                SM::HORadialIntegral(K-1,sps[mf].N,sps[mf].L,sps[mi].N,sps[mi].L,PhaseConvention);//now I multiply by Rfi
        return MagneticSPMultipole(MUM, sps, K, kappa, Rfi, mun, mup,GL_neutron,GL_proton);
    }
    
    
    
    
    /** \brief Single-particle Electric Multipole
     \f[
     \langle 2||{\cal M}(E\lambda)||1 \rangle = e \, i^{l_1-l_2}\,(-1)^{j_1-j_2+\lambda}
     \sqrt{\frac{(2\lambda+1)(2j_1+1)}{4\pi}}
     \langle j_1 \frac{1}{2}, \lambda 0 |j_2 \frac{1}{2} \rangle
     \langle 2|r^\lambda|1\rangle
     \f]
     Parity \f$ l_1+\lambda-l_2\f$ must be even.
     See also Bohr-Mottelson vol I, Eq. 3C-33, page 387.
     
     With this definition the reduced transition rate
     \f[
     B(E\lambda, I_1\rightarrow I_2)=\sum_{\mu, M_2} |\langle I_2 M_2|{\cal M}_{\lambda \mu}(E\lambda) | I_1 M_1\rangle |^2=
     \frac{\langle I_2||M_\lambda(E\lambda)||I_1\rangle}{2I_1+1}
     \f]
     The quadrupole operator is \f[Q=\sqrt{\frac{16\pi}{5}} {\cal M}_{\lambda \mu=0} \f]
     See also Bohr-Mottelson vol I, Eq. 3-30 , page 333.
     */
    template <class Radial>
    int ElectricSPMultipole(
                            MBOperator &MUM, ///<Multipole operator return
                            SingleParticleStates &sps, ///<Single particle states
                            int K, ///<Angular momentum
                            int kappa, ///<Magnetic projection
                            Radial &Rfi, ///<Radial Wave functions \f$ r^{K} \f$
                            double eneut, ///< neutron effective charge
                            double eprot ///<proton effective charge
    ){
        spsint pp[2]; //we first do pairs, must use constants
        double htmp;
        for (int mf=0;mf<sps.size();mf++) for (int mi=0;mi<sps.size();mi++) {
            //check parity and other simple stuff
            if  (sps[mf].Tz!=sps[mi].Tz) continue;
            if (((sps[mi].L+K-sps[mf].L)%2!=0)) continue;
            if ((sps[mi].Jz+2*kappa-sps[mf].Jz)!=0) continue;
            //reduced matrix element, page 387, volume I
            if (sps[mf].Tz==1) htmp=eneut; else htmp=eprot; //effective charge
            
            if (((((sps[mi].J-sps[mf].J)>>1)+K)%2)!=0) htmp=-htmp;
            //this was the factor of i^{li-lf+K} in front. is removed
//            if (((sps[mi].L - sps[mf].L)%4)!=0) htmp=-htmp;
            
            htmp*=sqrt((2.0*K+1.)*(sps[mi].J+1.)/(4.0*PI));
            htmp*=Rfi[sps[mf].level][sps[mi].level]; //now I multiply by Rfi
            htmp*= tav::ClebschGordan(sps[mi].J, 1,
                                      (2*K), 0,
                                      sps[mf].J, 1 );
            //if (std::abs(htmp)<1E-10)  continue; //zero due to radial happens often
            if (htmp==0)  continue;
            //reduced matrix element is done
            pp[0]=mi; pp[1]=mf;
//            MUM[pp] += htmp;
            MUM[pp]+=ReducedToFull(sps[mf],sps[mi],K,htmp);
        } //end loop
        
        return 0;
    }
    
    ///\brief Electric Multipole operator in Harmonic Oscillator wf
    /*
     (1/5/2016)--> Added a phase convention bool variable which is true by default but when set to true assumes wavefunctions that are positive at infinity (extra factor of \f$ (-1)^{n_1 + n_2} \f$

     */
    int ElectricSPMultipole(
                            MBOperator &MUM, ///<Multipole operator return
                            SingleParticleStates &sps, ///<Single particle states
                            int K, ///<Angular momentum
                            int kappa, ///<Magnetic projection
                            double eneut, ///< neutron effective charge
                            double eprot, ///<proton effective charge
                            bool PhaseConvention=true){
        
        // Determine the number of levels
        int levels=0;
        for (int mf=0;mf<sps.size();mf++) if (levels<sps[mf].level) levels=sps[mf].level;
        levels++;
        // create matrix
        tav::Tensor<2,double> Rfi(levels,levels);
        for (int mf=0;mf<sps.size();mf++)
            for (int mi=0;mi<sps.size();mi++) {
                //cerr<<"Multipole is started "<<levels<<" \n "<<sps[mf].N<<endl;
                Rfi[sps[mf].level][sps[mi].level]=
                SM::HORadialIntegral(K,sps[mf].N,sps[mf].L,sps[mi].N,sps[mi].L,PhaseConvention);//now I multiply by Rfi
            }
        return ElectricSPMultipole(MUM, sps, K, kappa, Rfi, eneut, eprot);
    }
    
    
    /*
     check for <1/2||s||1/2> was OK
     */
    /**
     \brief Spin-dependent operator
     \f[
     \langle 2||(Y_\kappa \sigma)_{(\kappa 1)\lambda}||1\rangle = i^{l_1-l_2} \sqrt{\frac{2j_1+1}{4\pi}} \left \{
     \begin{array}{l}
     (-1)^{j_1+\lambda-j_2}\sqrt{\lambda+1}\,\langle j_1 \frac{1}{2}, \lambda 0 |j_2 \frac{1}{2} \rangle - (-1)^{j_2-l_2-1/2}\sqrt{\lambda} \langle j_1 -\frac{1}{2}, \lambda 1 |j_2 \frac{1}{2} \rangle \quad \lambda=k-1 \\
     (-1)^{j_2-l_2-1/2} \sqrt{2\lambda+1}\,\langle j_1 -\frac{1}{2}, \lambda 1 |j_2 \frac{1}{2} \rangle \quad \lambda=k \\
     (-1)^{j_1+\lambda+j_2}\sqrt{\lambda}\langle j_1 \frac{1}{2}, \lambda 0 |j_2 \frac{1}{2} \rangle - (-1)^{j_2-l_2-1/2}\sqrt{\lambda+1} \,\langle j_1 -\frac{1}{2}, \lambda 1 |j_2 \frac{1}{2} \rangle \quad \lambda=k+1
     \end{array}
     \right .
     \f]
     See also Bohr-Mottelson I 3A-17 (p364)
     */
    double ROperator_YS(
                        QN &J2, ///< Final state quantum numbers: state \e 2
                        QN &J1, ///< Initial quantum numbers of state \e 1
                        int q, /// Angular momentum of Spherical Harmonics
                        int K /// Total momentum of the operator \f$ \lambda \f$
    ){
        
        if ((J1.P-J2.P+q)&1) return 0.0; //l1-l2+q is odd by parity requirement
        double htmp,htmp1, htmp2;
        //start with the common part
        htmp=sqrt((J1.J+1.)/(PI))/2.0; //common part includes 4pi and NO 1/2 for spin, use sigma
        //htmp*=PARITY(J1.L-J2.L+K); //this is a test pairity line not to be used
        switch (K-q) {
            case -1:
                htmp1=PARITY((((J1.J-J2.J)>>1)+K));
                htmp1*=tav::ClebschGordan(J1.J, 1,
                                          (2*K), 0,
                                          J2.J, 1 );
                htmp1*=sqrt(2.0*K+1.);
                //end of the first term
                htmp2= tav::ClebschGordan(J1.J, -1,
                                          (2*K), 2,
                                          J2.J, 1 );
                htmp2*=PARITY(((J2.J-1)>>1)-J2.L)*sqrt(double(K));
                htmp*=(htmp1-htmp2);
                return htmp;
            case 0:
                htmp*=PARITY(((J2.J-1)>>1)-J2.L);   //j_2-l_2-1/2
                htmp*=tav::ClebschGordan(J1.J, -1,
                                         (2*K), 2,
                                         J2.J, 1 );
                htmp*=sqrt(2.0*K+1.);
                return htmp;
            case 1: //this is GT operator
                htmp1=PARITY((((J1.J+J2.J)>>1)+K));
                htmp1*=tav::ClebschGordan(J1.J, 1,
                                          (2*K), 0,
                                          J2.J, 1 );
                htmp1*=sqrt(double(K));
                //end of the first term
                htmp2= tav::ClebschGordan(J1.J, -1,
                                          (2*K), 2,
                                          J2.J, 1 );
                htmp2*=PARITY(((J2.J-1)>>1)-J2.L)*sqrt(double(K)+1.0);
                htmp*=(htmp1-htmp2);
                return htmp;
        }
        return 0.0;
    }
    
    
    /** \brief Axial vector multipole moment
     \f[\langle 2||{\cal M}_A(\kappa, \lambda \mu)||1\rangle \,,
     \quad
     {\cal M}_A(\kappa, \lambda \mu) =\sum_p t^-(p) r_p^\kappa \left[Y_\kappa(p) \sigma(p) \right ]_{\lambda \mu}
     \f]
     
     See Bohr-Mottelson Eq. 3D-30, page 407.
     For evaluation use SM::ROperator_YS() function.
     
     The \f$ \kappa=0 \,\,\, \lambda=1 \f$ is Gamow-Teller transition \f$ {\cal M}_A(0, 1 \mu) =\frac{1}{\sqrt{4\pi}}\,\sum_p t^-(p)\sigma\mu(p) \f$
     
     Note,  the operator has no coupling constant \f$ g_a=(-1.23\pm 0.01)g_v \f$, for transition rate see BM 3D-38, page 411.
     \return Transition without \f$ g_a/\sqrt{4\pi} \f$ factor
     
     */
    template <class Radial>
    int BetaJAMultipole(
                        MBOperator &MUM, ///< output multipole matrix
                        SingleParticleStates &sps, ///< Single particle states
                        int K, ///< \f$ \lambda \f$ multipolarity of transition operator
                        int kappa, ///< magnetic projection
                        Radial &Rfi, ///< radial wave functions \f$ \langle 2|r^q|1\rangle \f$
                        int q  ///< Moment of Spherical Harmonics \f$ \kappa \f$
    ){
        /*
         
         NOTE OPERATOR HAS NO g_a/\sqrt{4\pi}
         
         */
        spsint pp[2];
        double htmp;
        for (int mf=0;mf<sps.size();mf++) {
            for (int mi=0;mi<sps.size();mi++) {
                //check parity and other simple stuff
                if  (sps[mf].Tz!=sps[mi].Tz+2) continue;  //T- operator
                if (((sps[mi].L-sps[mf].L+K)%2!=0)) continue; //(-1)^K is parity
                if ((sps[mi].Jz+2*kappa-sps[mf].Jz)!=0) continue;
                //reduced matrix element, page 387, volume I
                htmp=ROperator_YS(sps[mf], sps[mi], q, K);
                
                htmp*=Rfi[sps[mf].level][sps[mi].level]; //now I multiply by Rfi
                
                //cout<<htmp<<endl;
                
                
                //reduced matrix element is done
                //full matrix element, eq 1A-60, page 82
                /*
                 htmp*= ClebschGordan(sps[mi].J/2., sps[mi].Jz/2.,
                 double(K), double(kappa),
                 sps[mf].J/2., sps[mf].Jz/2.
                 );
                 
                 htmp/=sqrt(double(sps[mf].J+1));
                 */
                //if (std::abs(htmp)<1E-10)  continue; //zero due to radial happens often
                if (htmp==0)  continue;
                htmp=ReducedToFull(sps[mf],sps[mi],K,htmp);
                
                //also need t- operator
                htmp*=sqrt((sps[mi].T+sps[mi].Tz)*(sps[mi].T-sps[mi].Tz+2.))/2.0;
                
                //   if ((mi==mf)&&(sps[mi].J==sps[mi].Jz)) cerr<<sps[mi].Tz<<" "<<sps[mi].J/2.0<<" "<<htmp*sqrt(4.0*PI/3.0)<<endl; //for M1 these are sp moments
                pp[0]=mi; pp[1]=mf;
                MUM[pp]+=htmp;
            } }//end loop
        
        return 0;
    }
    
    /** \brief Beta multipole with Harmonic oscillator radial functions
     We produce combination of operators \f$ \tau^+ + \tau^- \f$
     they never both work at the same time so no-harm done
     Equivalent to version with general radial wave-functions.
     \todo discuss absent \f$ g_a/\sqrt{4\pi} \f$ factor
     
     (1/5/2016)--> Added a phase convention bool variable which is true by default but when set to true assumes wavefunctions that are positive at infinity (extra factor of \f$ (-1)^{n_1 + n_2} \f$

     */
    
    int BetaJAMultipole(
                        MBOperator &MUM,
                        SingleParticleStates &sps,
                        int K,
                        int kappa,
                        int q,
                        bool PhaseConvention=true){
        
        double htmp;
        spsint pp[2];
        
        for (int mf=0;mf<sps.size();mf++) {
            for (int mi=0;mi<sps.size();mi++) {
                //check parity and other simple stuff
                if  (sps[mf].Tz==sps[mi].Tz) continue;  //T- or T+ operator
                if (((sps[mi].L-sps[mf].L+q)%2!=0)) continue; //(-1)^q is parity Ja
                if ((sps[mi].Jz+2*kappa-sps[mf].Jz)!=0) continue;
                //reduced matrix element, page 387, volume I
                htmp=ROperator_YS(sps[mf], sps[mi], q, K);
                htmp*=SM::HORadialIntegral(q,sps[mf].N,sps[mf].L,sps[mi].N,sps[mi].L,PhaseConvention);//now I multiply by Rfi
                //if (std::abs(htmp)<1E-10)  continue; //zero due to radial happens often
                if (htmp==0)  continue;
                //cout<<htmp<<endl;
                
                
                //reduced matrix element is done
                //full matrix element, eq 1A-60, page 82
                /*
                 htmp*= ClebschGordan(sps[mi].J/2., sps[mi].Jz/2.,
                 double(K), double(kappa),
                 sps[mf].J/2., sps[mf].Jz/2.
                 );
                 
                 htmp/=sqrt(double(sps[mf].J+1));
                 */
                htmp=ReducedToFull(sps[mf],sps[mi],K,htmp);
                //also need t- operator
                if  (sps[mf].Tz==sps[mi].Tz-2) htmp*=sqrt((sps[mi].T+sps[mi].Tz)*(sps[mi].T-sps[mi].Tz+2.))/2.0; //T-
                else htmp*=sqrt((sps[mi].T-sps[mi].Tz)*(sps[mi].T+sps[mi].Tz+2.))/2.0;
                //T+
                pp[0]=mi; pp[1]=mf;
                MUM[pp]+=htmp;
            } }//end loop
        
        return 0;
    }
    
    /** \brief Beta plus  or minus multipole in HO basis 
     (1/5/2016)--> Added a phase convention bool variable which is true by default but when set to true assumes wavefunctions that are positive at infinity (extra factor of \f$ (-1)^{n_1 + n_2} \f$

     */
    int BetaJAMultipole(
                        MBOperator &MUM,
                        SingleParticleStates &sps,
                        int K,
                        int kappa,
                        int q,
                        int flag, ///< plus or minus flag=1 T+
                        bool PhaseConvention=true
    ){
//        cerr<<"beta (+/-) "<<flag<<" Phase Convention "<<PhaseConvention<<endl;
        double htmp;
        spsint pp[2];
        
        for (int mf=0;mf<sps.size();mf++) {
            for (int mi=0;mi<sps.size();mi++) {
                //check parity and other simple stuff
                if (flag==1) {if  (sps[mf].Tz!=sps[mi].Tz+2) continue;  }//T- or T+ operator
                else {if  (sps[mf].Tz!=sps[mi].Tz-2) continue;  }
                if (((sps[mi].L-sps[mf].L+q)%2!=0)) continue; //(-1)^q is parity Ja
                if ((sps[mi].Jz+2*kappa-sps[mf].Jz)!=0) continue;
                //reduced matrix element, page 387, volume I
                htmp=ROperator_YS(sps[mf], sps[mi], q, K);
                htmp*=SM::HORadialIntegral(q,sps[mf].N,sps[mf].L,sps[mi].N,sps[mi].L,PhaseConvention);//now I multiply by Rfi
                
                //cout<<htmp<<endl;
                
                
                //reduced matrix element is done
                //full matrix element, eq 1A-60, page 82
                /*
                 htmp*= ClebschGordan(sps[mi].J/2., sps[mi].Jz/2.,
                 double(K), double(kappa),
                 sps[mf].J/2., sps[mf].Jz/2.
                 );
                 
                 htmp/=sqrt(double(sps[mf].J+1));
                 */
                //if (std::abs(htmp)<1E-10)  continue; //zero due to radial happens often
               // cerr<<sps[mf].J<<" "<<sps[mi].J<<" "<<htmp<<endl;
                if (htmp==0)  continue;
                htmp=ReducedToFull(sps[mf],sps[mi],K,htmp);
                //also need t- operator below only "else" works
                if  (sps[mf].Tz==sps[mi].Tz-2) htmp*=sqrt((sps[mi].T+sps[mi].Tz)*(sps[mi].T-sps[mi].Tz+2.))/2.0; //T-
                else htmp*=sqrt((sps[mi].T-sps[mi].Tz)*(sps[mi].T+sps[mi].Tz+2.))/2.0; 
                //T+
                pp[0]=mi; pp[1]=mf;
                MUM[pp]+=htmp;
            } }//end loop
        
        return 0;
    }
    
    /**! \brief Make trivial multipole out of two single-particle states
    Note if need with clebsch then myltiply by sqrt(2K+1)
    */
    //      int Multipole original name
      int Multipole(
                            MBOperator &MUM, ///<Multipole operator return
                            SingleParticleStates &sps, ///<Single particle states
                            int K, ///<Angular momentum
                            int kappa, ///<Magnetic projection
                            int lev_c, ///< level id for creation operator
                            spin Tz_c, //creation operator
                            int lev_a, ///< level id for annihilation operator
                            spin Tz_a ///< annihilaiton operator isospin
    ){
        spsint pp[2]; //we first do pairs, must use constants
        double htmp;
        for (int mf=0;mf<sps.size();mf++) for (int mi=0;mi<sps.size();mi++) {
            //check parity and other simple stuff
            if ((sps[mi].Jz+2*kappa-sps[mf].Jz)!=0) continue;
            if (sps[mi].level!=lev_a)  continue;
            if (sps[mf].level!=lev_c) continue;
            if (sps[mi].Tz!=Tz_a)  continue;
            if (sps[mf].Tz!=Tz_c) continue;
            //reduced matrix element, page 387, volume I
            htmp=1.0;
            //reduced matrix element is done
            pp[0]=mi; pp[1]=mf;
            MUM[pp]+=ReducedToFull(sps[mf],sps[mi],K,htmp);
            } //end loop
        
        return 0;
    }

/*
This operator is   B dagger
\begin{equation}
 B^\dagger_{k \mu}(12)=\sum_{m_1 m_2} (-)^{j_1-m_1}  \begin{pmatrix}
    j_1 & k & j_2 \\
    -m_1 & \mu & m_2
  \end{pmatrix}
  a^\dagger_{j_1 m_1} a_{j_2 m_2}
 \end{equation}
 Equivalently we can express this using CGC
 \begin{equation}
 B^\dagger_{k \mu}(12)=\sum_{m_1 m_2} \frac{(-)^{2k}}{\sqrt{2j_1+1}} C^{j_1 m_1}_{j_2 m_2\, k \mu}  a^\dagger_{j_1 m_1} a_{j_2 m_2}
 \end{equation}
\begin{equation}
 B^\dagger_{k \mu}(12)=\sum_{m_1 m_2} \frac{(-)^{2j_1}}{\sqrt{2k+1}} (-)^{j_2-m_2} C^{k \mu}_{j_1 m_1\, j_2 m_2}  a^\dagger_{j_1 m_1} a_{j_2 -m_2}
 \end{equation}
(in the book  ${\cal M}^\dagger_{k\mu}= \sqrt{2k+1} B^\dagger_{k \mu},$
*/
//this is uneccessary repeates above
//    int ParticleHoleOperator(
//                            MBOperator &MUM, ///<Multipole operator return
//                            SingleParticleStates &sps, ///<Single particle states
//                            int K, ///<Angular momentum
//                            int kappa, ///<Magnetic projection
//                            int l1, ///<level index 1
//                            int l2, ///<level index 2
//                            double eneut, ///< neutron effective charge
//                            double eprot ///<proton effective charge
//    ){
//        spsint pp[2]; //we first do pairs, must use constants
//        double htmp;
//        for (int mf=0;mf<sps.size();mf++) for (int mi=0;mi<sps.size();mi++) {
//            //check parity and other simple stuff
//            if (sps[mf].level!=l1) continue;
//             if (sps[mi].level!=l2) continue;
//            if  (sps[mf].Tz!=sps[mi].Tz) continue;
//           // if (((sps[mi].L+K-sps[mf].L)%2!=0)) continue;
//            if ((sps[mi].Jz+2*kappa-sps[mf].Jz)!=0) continue;
//            //reduced matrix element, page 387, volume I
//            if (sps[mf].Tz==1) htmp=eneut; else htmp=eprot; //effective charge
//            
//            //if (std::abs(htmp)<1E-10)  continue; //zero due to radial happens often
//            if (htmp==0)  continue;
//            //reduced matrix element is done
//            pp[0]=mi; pp[1]=mf;
////            MUM[pp] += htmp;
////            MUM[pp]+=htmp*tav::ThreeJSymbol(JI.J, JI.Jz,
////                                            2*Lambda, (JF.Jz-JI.Jz),
////                                            JF.J, JF.Jz
////                                            )/sqrt(double(JF.J+1));
//            
//           // MUM[pp]+=SM::ReducedToFull(sps[mf],sps[mi],K,htmp);
//            MUM[pp]+=htmp*tav::ClebschGordan(sps[mi].J, sps[mi].Jz,
//                                                2*K, (sps[mf].Jz-sps[mi].Jz),
//                                                sps[mf].J, sps[mf].Jz
//                                                )/sqrt(double(sps[mf].J+1));
//        } //end loop
//        
//        return 0;
//    }

    /**! \brief Make trivial pair two single-particle states
    Note this we think of creation
     \begin{equation}
 A^\dagger_{k \mu}(12)=\sum_{m_1 m_2} (-)^{j_1-m_1+j_2-m_2}  
 \begin{pmatrix}
    j_1 & k & j_2 \\
    -m_1 & \mu & m_2
  \end{pmatrix}
  a^\dagger_{j_1 m_1} a^\dagger_{j_2 -m_2},
 \end{equation}
\begin{equation}
A^\dagger_{k \mu}(12)=\sum_{m_1 m_2} \frac{(-)^{2k}}{\sqrt{2k+1}} C^{k \mu}_{j_1 m_1\, j_2 m_2}  a^\dagger_{j_1 m_1} a^\dagger_{j_2 m_2}
 \end{equation}
 Therefore for a singe level also have $P^\dagger_{LM}=\sqrt{\frac{2L+1}{2}}A^\dagger_{LM}.$ 
 
    */
      int PairOperator(
                            MBOperator &MUM, ///<Pair operator return
                            SingleParticleStates &sps, ///<Single particle states
                            int K, ///<Angular momentum
                            int kappa, ///<Magnetic projection
                            int lev_1, ///< level id for left creation operator
                            spin Tz_1, //isospin of left
                            int lev_2, ///< level id for right creation operator
                            spin Tz_2 ///< isospin of right
    ){
        spsint pp[2]; //we first do pairs, must use constants
        double htmp;
        for (int m1=0;m1<sps.size();m1++) for (int m2=0;m2<sps.size();m2++) {
            //check parity and other simple stuff
            if (m1==m2) continue;
            if ((sps[m2].Jz-2*kappa+sps[m1].Jz)!=0) continue;
            if (sps[m2].level!=lev_2)  continue;
            if (sps[m1].level!=lev_1) continue;
            if (sps[m2].Tz!=Tz_2)  continue;
            if (sps[m1].Tz!=Tz_1) continue;
            htmp=1.0;
            //reduced matrix element is done
            
            if (m1<m2) {
            pp[0]=m1; pp[1]=m2;
            MUM[pp]+=
            tav::ClebschGordan(
            sps[m1].J, sps[m1].Jz,
            sps[m2].J, sps[m2].Jz,
            2*K, 2*kappa
            );}
            else {
            pp[0]=m2; pp[1]=m1;
            MUM[pp]-=
            tav::ClebschGordan(
            sps[m1].J, sps[m1].Jz,
            sps[m2].J, sps[m2].Jz,
            2*K, 2*kappa
            );
            }

            
           // ReducedToFull(sps[mf],sps[mi],K,htmp);
            } //end loop
        
        return 0;
    }



        /*
     Creation/annihilarion operator for CM quanta
     (1/5/2016)--> Added a PhaseConvention bool variable which is false by default
     */
    int CMCreationOperator(
                           MBOperator &MUM, ///<Multipole operator return
                           SingleParticleStates &sps, ///<Single particle states
                           int otype, //type +1 for creation -1 for annihilation
                           int kappa, ///<Magnetic projection
                           double eneut, ///< neutron effective charge
                           double eprot, ///<proton effective charge
                            bool PhaseConvention=true
    ){
        int K=1;//multipolarity
        spsint pp[2]; //we first do pairs, must use constants
        double htmp;
        for (int mf=0;mf<sps.size();mf++) for (int mi=0;mi<sps.size();mi++) {
            //check parity and other simple stuff
            if  (sps[mf].Tz!=sps[mi].Tz) continue;
            if (((sps[mi].L+K-sps[mf].L)%2!=0)) continue;
            if ((sps[mi].Jz+2*kappa-sps[mf].Jz)!=0) continue;
            if (2*sps[mf].N+sps[mf].L!=2*sps[mi].N+sps[mi].L+otype) continue;
            //reduced matrix element, page 387, volume I
            if (sps[mf].Tz==1) htmp=eneut; else htmp=eprot; //effective charge
            if (((((sps[mi].J-sps[mf].J)>>1)+K)%2)!=0) htmp=-htmp;
            //This was the i^(li+K-lf) part. commented out now.
//            if (((sps[mi].L-sps[mf].L)%4)!=0) htmp=-htmp;
            htmp*=sqrt((2.0*K+1.)*(sps[mi].J+1.)/(4.0*PI));
            htmp*=SM::HORadialIntegral(K,sps[mf].N,sps[mf].L,sps[mi].N,sps[mi].L,PhaseConvention);
            //Rfi[sps[mf].level][sps[mi].level]; //now I multiply by Rfi
            htmp*= tav::ClebschGordan(sps[mi].J, 1,
                                      (2*K), 0,
                                      sps[mf].J, 1 );
            //if (std::abs(htmp)<1E-10)  continue; //zero due to radial happens often
            if (htmp==0)  continue;
            //reduced matrix element is done
            pp[0]=mi; pp[1]=mf;
            MUM[pp]+=ReducedToFull(sps[mf],sps[mi],K,htmp);
        } //end loop
        
        return 0;
    }
    //These are not used.
    //They where added in 2016-17 when checking the phases of transition operators.

//    double YLRMEj(int L, int lf, int jf, int li, int ji)//j,jp are half integer so these are double values
//    {
//        if ((lf+L+li)&1) return 0;
//        double rtrn = std::sqrt((jf+1) * (2*L+1) * (ji+1)/ 4. /M_PI);
//        rtrn *= mav::ThreeJSymbol(jf/2.,-0.5, L,0, ji/2.,0.5);
//        if ((((jf - 1)>>1) +L )   &1) rtrn = -rtrn;//j-1/2 phase
//        if ((((1-jf)>>1) +lf )&1) rtrn = -rtrn;//phase for difference in order of coupling. (l+s=j)
//        if ((((1-ji)>>1) +li )&1) rtrn = -rtrn;//same for initial state.
//        //BM phases
////        if (((ji-jf)>>1 + L)&1) rtrn = -rtrn;
////        if (((li-lf+L)>>2)&1) rtrn = -rtrn;
//        return rtrn;
//    }
//    double YLRME(int L,int lf, int li)
//    {
//        double rtrn = mav::ThreeJSymbol(lf,0.,L,0.,li,0);
//        rtrn *= std::sqrt((2*L+1) * (2*lf+1) * (2*li+1)/4./M_PI);
//        if (lf&1) rtrn = -rtrn;
//        return rtrn;
//    }
//    template <class Radial>
//    int ElectricSPMultipole2(
//                             MBOperator &MUM, ///<Multipole operator return
//                             SingleParticleStates &sps, ///<Single particle states
//                             int K, ///<Angular momentum
//                             int kappa, ///<Magnetic projection
//                             Radial &Rfi, ///<Radial Wave functions \f$ r^{K} \f$
//                             double eneut, ///< neutron effective charge
//                             double eprot ///<proton effective charge
//    ){
//        spsint pp[2]; //we first do pairs, must use constants
//        double htmp;
//        for (int mf=0;mf<sps.size();mf++)
//        {
//            for (int mi=0;mi<sps.size();mi++)
//            {
//                if (sps[mf].Tz!=sps[mi].Tz) continue;
//                if (sps[mi].Jz+2*kappa-sps[mf].Jz) continue;
//                if ((sps[mf].L+K+sps[mi].L)&1) continue;//parity. Included in RME of Y_L
//                if (sps[mi].J + (K<<1) < sps[mf].J) continue;//triangle rule (+)
//                if (sps[mi].J - (K<<1) > sps[mf].J) continue;//triangle rule (-)
//                
//                //now actual calculation.
//                htmp = (sps[mf].Tz==1 ? eneut : eprot); //effective charge
//                htmp *= Rfi[sps[mf].level][sps[mi].level]; //radial part. Includes case of (-)^n phase convention
//                if (std::abs(htmp) < 1.e-14) continue;
//                htmp *= YLRMEj(K,sps[mf].L,sps[mf].J,sps[mi].L,sps[mi].J);
//                if (std::abs(htmp) < 1.e-14) continue;
//
//                //reduced to full
//                htmp *= mav::ThreeJSymbol(sps[mf].J/2., -sps[mf].Jz/2.,
//                                          K,kappa,
//                                          sps[mi].J/2., sps[mi].Jz/2.);
////                if (((sps[mf].Jz+sps[mf].J)>>1+sps[mi].J)&1) htmp = -htmp;
//                if (((sps[mf].J - sps[mf].Jz)>>1)&1) htmp = -htmp;
//                if (std::abs(htmp) < 1.e-14) continue;
//                
////                if ((sps[mf].L + (sps[mf].L + sps[mi].L)/2)&1) htmp = -htmp;
////                if ((sps[mi].L - sps[mf].L +K)%4 != 0) htmp = -htmp;
////                if (((sps[mi].J-1)>>1)&1) htmp = -htmp;
//
//                //and store.
//                pp[0]=mi; pp[1]=mf;
//                MUM[pp] += htmp;
//            }
//        }
//        return 0;
//    }
//    int ElectricSPMultipole2(
//                            MBOperator &MUM, ///<Multipole operator return
//                            SingleParticleStates &sps, ///<Single particle states
//                            int K, ///<Angular momentum
//                            int kappa, ///<Magnetic projection
//                            double eneut, ///< neutron effective charge
//                            double eprot, ///<proton effective charge
//                            bool PhaseConvention=true){
//        
//        // Determine the number of levels
//        int levels=0;
//        for (int mf=0;mf<sps.size();mf++) if (levels<sps[mf].level) levels=sps[mf].level;
//        levels++;
//        // create matrix
//        tav::Tensor<2,double> Rfi(levels,levels);
//        for (int mf=0;mf<sps.size();mf++)
//            for (int mi=0;mi<sps.size();mi++) {
//                //cerr<<"Multipole is started "<<levels<<" \n "<<sps[mf].N<<endl;
//                Rfi[sps[mf].level][sps[mi].level]=
//                SM::HORadialIntegral(K,sps[mf].N,sps[mf].L,sps[mi].N,sps[mi].L,PhaseConvention);//now I multiply by Rfi
//            }
//        return ElectricSPMultipole2(MUM, sps, K, kappa, Rfi, eneut, eprot);
//    }
    
} //end namespace

/**
 \page angular Angular Momentum Recoupling
 See also A.R. Edmonds, Angular Momentum in Quantum Mechanics. 
 \par Reduced Matrix Element
 In the angular momentum the operator behaves as adding j'=\lambda+j we define
 a reduced matrix element as
 \f[
 \frac{\langle j'm'|Q_{\lambda\mu}^{\dagger}|jm\rangle}{\langle j'||Q_{\lambda}^{\dagger}||j\rangle}=(-)^{j'-m'}\left(\begin{array}{ccc}
 j' & \lambda & j\\
 -m' & \mu & m\end{array}\right)=\frac{(-)^{j-m}}{\sqrt{2\lambda+1}}C_{j'm'j-m}^{\lambda\mu}=\frac{(-)^{\lambda-j+j'}}{\sqrt{2j'+1}}C_{\lambda\mu jm}^{j'm'}=\frac{(-)^{2\lambda}}{\sqrt{2j'+1}}C_{jm\lambda\mu}^{j'm'}
 \f] 
 see [Edmonds 5.4.1]
 \par  Wigner 3j-symbols
 
 The relation between Clebsch-Gordan coefficients and 3j-symbols: 
 \f[ 
 C_{j_{1}m_{1}\; j_{2}m_{2}}^{j_{3}m_{3}}=(-)^{j_{1}-j_{2}+m_{3}}\sqrt{2j_{3}+1}\,\left(\begin{array}{ccc}
 j_{1} & j_{2} & j_{3}\\
 m_{1} & m_{2} & -m_{3}\end{array}\right).
 \f]
 
 The simple cases correspond to angular momentum zero and one: 
 \f[ 
 \left(\begin{array}{ccc}
 j & j & 0\\
 m & -m' & 0\end{array}\right)=C_{jm\, jm'}^{00}=\delta_{m',-m}(-)^{j-m}\frac{1}{(2j+1)^{1/2}}, \left(\begin{array}{ccc}
 j & j & 1\\
 m & -m & 0\end{array}\right)=(-)^{j-m}\frac{M}{\sqrt{(2j+1)j(j+1)}}.
 \f]
 Orthogonality properties 
 \f[
 \sum_{j_{3}m_{3}}(2j_{3}+1)\left(\begin{array}{ccc}
 j_{1} & j_{2} & j_{3}\\
 m_{1} & m_{2} & m_{3}\end{array}\right)\left(\begin{array}{ccc}
 j_{1} & j_{2} & j_{3}\\
 m'_{1} & m'_{2} & m_{3}\end{array}\right)=\delta_{m_{1}m'_{1}}\delta_{m_{2}m'_{2}} \sum_{m_{1}m_{2}}\left(\begin{array}{ccc}
 j_{1} & j_{2} & j_{3}\\
 m_{1} & m_{2} & m_{3}\end{array}\right)\left(\begin{array}{ccc}
 j_{1} & j_{2} & j'_{3}\\
 m_{1} & m_{2} & m'_{3}\end{array}\right)=\frac{\delta_{j_{3}j'_{3}}\delta_{m_{3}m'_{3}}}{2j_{3}+1}\,\delta(j_{1}j_{2}j_{3}) 
 \f]
 
 \par Wigner 6j-symbols
 
 Simple cases: 
 \f[
 \left\{ \begin{array}{ccc}
 a & b & c\\
 0 & c' & b'\end{array}\right\} =\delta_{b'b}\delta_{c'c}\frac{(-)^{a+b+c}}{\sqrt{(2b+1)(2c+1)}}. \left\{ \begin{array}{ccc}
 a & b & c\\
 \frac{1}{2} & c-\frac{1}{2} & b\pm\frac{1}{2}\end{array}\right\} =(-)^{a+b+c}\sqrt{\frac{(\pm a\mp b+c)(a+b\mp c+1)}{(2b+1)(2c+1)2c(2b+1/2\pm1/2)}}.
 \f]
 The recoupling property:
 \f[
 \sum_{M'}(-)^{J'+M'}\left(\begin{array}{ccc}
 j_{1} & j_{2} & J'\\
 m_{1} & m_{2} & M'\end{array}\right)\left(\begin{array}{ccc}
 j_{3} & j_{4} & J'\\
 m_{3} & m_{4} & -M'\end{array}\right)= \sum_{JM}(-)^{2j_{4}+J+M}(2J+1)\left\{ \begin{array}{ccc}
 j_{1} & j_{2} & J'\\
 j_{3} & j_{4} & J\end{array}\right\} \left(\begin{array}{ccc}
 j_{3} & j_{2} & J\\
 m_{3} & m_{2} & M\end{array}\right)\left(\begin{array}{ccc}
 j_{1} & j_{4} & J\\
 m_{1} & m_{4} & -M\end{array}\right). 
 \f]
 Another form of the recoupling identity:
 \f[
 \sum_{\mu_{1}\mu_{2}\mu_{3}}(-)^{l_{1}+l_{2}+l_{3}+\mu_{1}+\mu_{2}+\mu_{3}}\left(\begin{array}{ccc}
 j_{1} & l_{2} & l_{3}\\
 m_{1} & \mu_{2} & -\mu_{3}\end{array}\right)\left(\begin{array}{ccc}
 l_{1} & j_{2} & l_{3}\\
 -\mu_{1} & m_{2} & \mu_{3}\end{array}\right)\left(\begin{array}{ccc}
 l_{1} & l_{2} & j_{3}\\
 \mu_{1} & -\mu_{2} & m_{3}\end{array}\right) =\left(\begin{array}{ccc}
 j_{1} & j_{2} & j_{3}\\
 m_{1} & m_{2} & m_{3}\end{array}\right)\left\{ \begin{array}{ccc}
 j_{1} & j_{2} & j_{3}\\
 l_{1} & l_{2} & l_{3}\end{array}\right\} .
 \f]
 The above recoupling can be expressed in the form of Clebsh-Gordan coefficients
 \f[
 \sum_{m_{12}}C_{j_{1}m_{1}j_{2}m_{2}}^{j_{12}m_{12}}C_{j_{12}m_{12}j_{3}m_{3}}^{jm}=\sum_{j_{23}m_{23}}(-)^{j_{1}+j_{2}+j_{3}+j}\sqrt{(2j_{12}+1)(2j_{23}+1)}\left\{ \begin{array}{ccc}
 j_{1} & j_{2} & j_{12}\\
 j_{3} & j & j_{23}\end{array}\right\} \times C_{j_{2}m_{2}j_{3}m_{3}}^{j_{23}m_{23}}C_{j_{1}m_{1}j_{23}m_{23}}^{jm}
 \f]
 The same recoupling for the particle-hole representation
 \f[
 \sum_{\mu}C_{j_{1}m_{1}\lambda\mu}^{j_{2}m_{2}}C_{j_{3}m_{3}\lambda\mu}^{j_{4}m_{4}}=\sum_{LM}(-)^{j_{1}+j_{2}+j_{3}+j_{4}}\sqrt{(2j_{2}+1)(2j_{4}+1)}\left\{ \begin{array}{ccc}
 j_{1} & j_{2} & \lambda\\
 j_{3} & j_{4} & L\end{array}\right\} \times C_{j_{3}m_{3}j_{2}m_{2}}^{LM}C_{j_{1}m_{1}j_{4}m_{4}}^{LM}
 \f]
 
 The orthogonality relations for the 6j-symbols are, 
 \f[
 (2l+1)\sum_{K}(2K+1)\left\{ \begin{array}{ccc}
 a & b & l\\
 c & d & K\end{array}\right\} \left\{ \begin{array}{ccc}
 c & b & K\\
 a & d & l'\end{array}\right\} =\delta_{ll'}, and \sum_{K}(-)^{K+l+l'}(2K+1)\left\{ \begin{array}{ccc}
 a & b & l\\
 c & d & K\end{array}\right\} \left\{ \begin{array}{ccc}
 b & c & K\\
 a & d & l'\end{array}\right\} =\left\{ \begin{array}{ccc}
 b & a & l\\
 c & d & l'\end{array}\right\} .
 \f]
 */


#endif 



