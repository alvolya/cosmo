# Example: study of 33Al with fsu Interaction

## 1) Run cosmoxml for System Description
**Command:** `cosmoxml -nn 33Al`  
**User Input:**  
- Select model space: `spsdpf`  
- Select Hamiltonian: `fsu9`  
- Enter spin projection Jz[1/2]: (press Enter)  
- Would you like to shift Hamiltonian matrix by J+J- term(y/n)[n]: (press Enter)  
- Enter parity (+ or -): `+`  
- Create hw rejection (use with Xsysmbs)(y/n)[n]: `y`  
- Maximum number of oscillator quanta above minimum: [0]  (press Enter) 
- Minimum number of oscillator quanta above minimum: [0]  (press Enter)
- Create Particle-Hole Rejection (use with XSHLMBSQRR)(y/n)[n]: (press Enter)  
- Create rejection (to be processed by Xsysmbs)(y/n)[n]: (press Enter)  
- Name the system as `[33Al_fsu9_M1T7+]`: (press Enter)  
- Name output xml file as `[33Al_fsu9+.xml]`: (press Enter)  

**Comment:** 
> This creates 0hw model space with positive parity states.

## Repeat for 1hw States with Negative Parity
**Comment:**  
> Similar steps for negative parity but different choices for parity and oscillator quanta (min 1 and max 1).  
> Comment: Additional systems for 2hw and 3hw states may be obtained similarly but require different system names and XML file names.


## 2) Create Basis for Both Parity States
**Commands:** 
- `XsysmbsOMP 33Al_fsu9+.xml`
- `XsysmbsOMP 33Al_fsu9-.xml`

**Program Output:**
- Number of mbs 119 for positive parity
- Number of mbs 23842 for negative parity

## 3) Diagonalize Both Systems
**Commands:** 
- `DEXHVJTcsb 33Al_fsu9+.xml`
- `DEXHVJTcsb 33Al_fsu9-.xml`

**Program Output for Positive Parity:**
- "There are 119, how many do you want to process?" User Input: `10` 
**Comment:** 
> because system is small full diagonalization was used and all 119 states are available (see -l flag). 

## 4) Examine the Resulting Spectrum
**Command:** `Xlevels 33Al_fsu9+.xml 33Al_fsu9-.xml`  
**Program Output:**
- 5/2+(1) T=7/2 0 -303.079
- 1/2-(1) T=7/2 1.947 -301.132
- 7/2-(1) T=7/2 3.151 -299.928
- 5/2-(1) T=7/2 3.205 -299.874
- 3/2-(1) T=7/2 3.262 -299.817

**Comment:** 
> Note the state at 3.967 MeV is the last state of negative parity and was not well resolved. Additional iterations in Davidson may be required.

## 5) Analyze E2 Transitions Between Positive Parity States
**Command:** `XSHLEMB 33Al_fsu9+.xml -u -w`  
**User Input:**
- Operator type: `E` for electric
- Multipolarity: `2` for E2
- Neutron effective charge: `0.5`   (enter neutron effective charge)
- Proton effective charge: `1.5`  
- File for radial overlap: `HO`  (enter HO for harmonic oscillator or file name for matrix of other radial wave function overlaps) 

**Program Output:**
- Transition details..
- Output stored in XML files with details, The B(E2) is given as  B is in oscillator units, Bu is in fm and weisskopf is for Weisskopf units (see XSHLEMB --help for details). 
