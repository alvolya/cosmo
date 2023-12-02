# 24Mg Calculation with usdb Interaction

## 1) Create an XML File with System Description 
**Command:** `> cosmoxml`  
**User Input:**  
- Enter nucleus name: `24Mg`  
- Select model space: `sd`  
- Select Hamiltonian: `usdb`  
- Spin projection Jz[0]: (press Enter for default `[0]`)  
- Shift Hamiltonian matrix by J+J- term (y/n)[n]: (press Enter for default `n`)  
- Create hw rejection (use with Xsysmbs)(y/n)[n]: (press Enter for default `n`)  
- Create Particle-Hole Rejection (use with XSHLMBSQRR)(y/n)[n]: (press Enter for default `n`)  
- Create rejection (to be processed by Xsysmbs)(y/n)[n]: (press Enter for default `n`)  
- Name the system as [24Mg_usdb_M0T0+]: (press Enter for default)  
- Name output xml file as [24Mg_usdb+.xml]: (press Enter for default)

**Comment:**  
> An XML file with system description is created. More about system description can be found at [www.volya.net - XML Format](https://www.volya.net/index.php?id=xml-format).

## 2) Create Basis, Diagonalize the Hamiltonian, determine spins of states, and populate xml file
**Command:** `DEXHVJTcsb 24Mg_usdb+.xml`

## 3) Review States
**Command:** `Xlevels 24Mg_usdb+.xml`  
**User Input:**  
- Include this space (y/n)[y]: (press Enter for default `y`)  

**Program Output:**  
> 0+(1) T=0 0 -87.1044  
> 2+(1) T=0 1.5023 -85.6021  
> 2+(2) T=0 4.1161 -82.9883  
> 4+(1) T=0 4.3724 -82.732  
> 3+(1) T=0 5.075 -82.0294  
> 4+(2) T=0 5.8825 -81.2219  
> 0+(2) T=0 7.3382 -79.7662  
> 2+(3) T=0 7.4817 -79.6227  
> 5+(1) T=0 7.7968 -79.3076  
> 1+(1) T=0 7.8181 -79.2863
