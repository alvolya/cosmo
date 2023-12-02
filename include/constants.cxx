/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


#ifndef __CONSTANTS__CXX
#define __CONSTANTS__CXX
/*! \namespace phy
\brief physical constants and functions
 */
/*! \file constants.cxx
\brief Physical Constants
*/
namespace phy{
    //! \brief \f$\hbar c=197.32\f$ MeV fm            
  const double HBARC=197.327053;
   //! \brief Atomic mass u=931.5 MeV, defined via mass of carbon 
  const double AMU=931.494043;
     //! \brief Proton mass in MeV
  const double PMASS=938.272029;
     //! \brief Neutron mass in MeV
  const double NMASS=939.565360; 
  //! \brief Electron mass in MeV  
  const double EMASS=0.510998910;
     //! \brief speed of light in meters per second
  const double C=299792458; 
   //! \brief \f$ \alpha=e^2/(\hbar c) \approx 1/137 \f$ fine structure constant
  const double alpha=7.297352568E-3;  //1/137 use of e^2
  //!magnetic moments of proton, in units of nuclear magneton \f$ e\hbar/(2M_p c)\f$ remember \f$ g=2\mu \f$
  const double muproton=2.792847351; 
//! magnetic moment of a neutron in units of nuclear magnetons, see #muproton
  const double muneutron=-1.91304273; 
  const double pi=3.141592653589793238462643;
//PI is already used in define, but small pi is for constant...
  const double euler_constant=0.577215664901532860606512; 
  //! Boltzmann constant in MeV per degree of K 
  const double Kb=8.617343E-11; //MeV K-1
}
#endif
