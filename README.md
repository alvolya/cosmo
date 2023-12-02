# CoSMo - Shell Model Code

CoSMo is a comprehensive shell model code suite designed for nuclear shell model and configuration interaction calculations. The highly structural and templated nature of the CoSMo code allows for flexibility and ease in applications, including those to open quantum systems with non-Hermitian Hamiltonians, clustering, and time-dependent dynamics.

## Installation Notes

Before installing, ensure your g++ compiler is compatible. Modify `Makefile` and `build.sh` as necessary to match your system's configuration.

To install:
1. Run `build.sh`. This script will compile all files into the `./bin` directory.
2. During the first run, the Eigen library will be installed via `git clone https://gitlab.com/libeigen/eigen.git`.

## Usage

Before using CoSMo, source the appropriate environment variables with either `source cosmovar.sh` or `source cosmovar.zsh`.

To (re)compile an individual `file.cpp`, use the command `make file`.

## Running CoSMo

See EXAMPLE* files provided, additional information can be found at www.volya.net

1. **Create XML File**: 
   - Create an XML file that identifies your system. This can be done manually, by editing a copied template, or using `cosmoxml` or `xcosmo` which interact with the database of all interactions.
   - Command: `cosmoxml`

2. **Create Many-Body Basis**: 
   - Use `Xsysmbs myfile.xml` or `XsysmbsOMP` for the OpenMP version.

3. **Create Many-Body Hamiltonian** (Optional if using `DEXHVJTcsb`): 
   - Command: `XHH+JJ myfile.xml`

4. **Diagonalize Hamiltonian** (Optional if using `DEXHVJTcsb`): 
   - Use `davidson_file myfile.HH` for large matrices or `exactev` (`texactev` or `texactevEigen`) for full diagonalization.

5. **Determine Spin and Isospin** (Optional if using `DEXHVJTcsb`): 
   - Command: `XSHLJT myfile.xml`

6. **Compute EM Transitions**: 
   - Command: `XSHLEMB myfile.xml -o output.xml`

7. **Compute Spectroscopic Factors**: 
   - Command: `XSHLSF myfile1.xml -f myfile2.xml`

8. **Compute Occupation Numbers**: 
   - Command: `XSHLAO myfile.xml`

9. **Alternative Combined Command**: 
   - `DEXHVJTcsb myfile.xml` can replace steps 3, 4, and 5.

10. **Viewing and Editing XML Database**: 
    - Use `Xlevels myfile.xml` to see an organized list of levels.
    - XML files are regular text files and can be edited directly for modifications or to set up new jobs.

## Basic Codes

The files included in the CoSMo distribution are listed in `cosmo-main.txt`. All codes should accept command-line arguments and support `-h` or `--help` for help and usage instructions.

- `cosmoxml`: Generate an XML file for a system.
- `Xsysmbs`: Create many-body states (*.mbs files), with support for rejections.
- `XsysmbsOMP`: OpenMP version of `Xsysmbs`.
- `DEXHVJT`: Diagonalize Hamiltonian virtually and store states in an XML file.
- `DEXHVJTcsb`: An alternative method for diagonalizing the Hamiltonian.
  - Combines `XHH+JJ`, `XSHLJT`, and `davidson_vfile`.
- `XHH+JJ.cpp`: Create a Hamiltonian matrix HH.
- `davidson_vfile.cpp`: Diagonalize the Hamiltonian matrix.
- `XSHLJT.cpp`: Process eigenvectors, determine spin, and populate the XML file.
- `Xlevels`: Display states from the XML file.
- `XSHLEMB`: Study electromagnetic and beta transitions.
- `XSHLSF`: Study classic SF (no recoil).
- `XSHLAO`: View occupation numbers.


## License


This project is licensed under the [GNU General Public License v3 (GPL v3)](./LICENSE). For the full license text, see the `LICENSE` file.

For more details and acknowledgments of third-party components and libraries, please refer to the [`NOTICE`](./NOTICE) file.
