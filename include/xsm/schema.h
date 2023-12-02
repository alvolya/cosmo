/*! \file schema.h
\ingroup gp_sm
*/
#ifndef __SCHEMA__H__
#define __SCHEMA__H__
namespace xsm{
const char cosmo[]="cosmo";
const char levels[]="levels";
const char system[]="system";
const char file[]="file";
const char nucleus[]="nucleus"; ///< used in exp data xlevels
const char system_basisstates[]="ManyBodyStates";
const char system_matrix[]="matrix";
const char system_gmatrix[]="G-matrix";
const char system_eigenstate[]="state";
const char system_channel[]="channel";
const char system_eigenstate_wavefunction[]="wavefunction";
const char system_operator[]="operator";


const char hamiltonian[]="hamiltonian";
const char hamiltonian_comment[]="comment";
//const char hamiltonian_reference[]="reference";
const char hamiltonian_formulascale[]="formulascale";
//const char hamiltonian_file[]="file";
//const char hamiltonian_file_name[]="name";
const char JISP[]="jispfile";

const char valencespace[]="valencespace";
const char valencespace_orbital[]="orbital";
const char valencespace_core[]="core";
const char valencespace_outercore[]="outercore";





//attributes
//QN
const char J[]="J";
const char Jz[]="Jz";
const char T[]="T";
const char Tz[]="Tz";
const char P[]="P";
const char E[]="E";
const char Ex[]="Ex";
const char name[]="name";
const char Nmax[]="Nmax"; //max number of oscillator quanta
const char hwmax[]="hwmax"; //max number of oscillator excitation quanta
const char maxshell[]="maxshell"; //maximum number of quanta in sp state
//generally speaking Nmax and maxshell must be the same 

//valence
const char valence_A[]="Av";
const char valence_N[]="Nv"; //currently not used
const char valence_Z[]="Zv"; //currently not used
//nucleus
const char nucleus_A[]="A";
const char nucleus_N[]="N";
const char nucleus_Z[]="Z";

}
#endif //__SCHEMA__H__
