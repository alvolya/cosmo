/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


#ifndef __NUCLEARNAMES_CXX
#define __NUCLEARNAMES_CXX
#include <debug.h>
#include <cstring>
#include <cstdlib> //atoi
namespace phy{
char NuclearSymbol[119][4]={"n",  "H",  "He",  "Li",  "Be",   
"B",  "C",  "N",  "O",  "F",   
"Ne",  "Na",  "Mg",  "Al",  "Si",   
"P",  "S",  "Cl",  "Ar",  "K",   
"Ca",  "Sc",  "Ti",  "V",  "Cr",   
"Mn",  "Fe",  "Co",  "Ni",  "Cu",   
"Zn",  "Ga",  "Ge",  "As",  "Se",   
"Br",  "Kr",  "Rb",  "Sr",  "Y",   
"Zr",  "Nb",  "Mo",  "Tc",  "Ru",   
"Rh",  "Pd",  "Ag",  "Cd",  "In",   
"Sn",  "Sb",  "Te",  "I",  "Xe",   
"Cs",  "Ba",  "La",  "Ce",  "Pr",   
"Nd",  "Pm",  "Sm",  "Eu",  "Gd",   
"Tb",  "Dy",  "Ho",  "Er",  "Tm",   
"Yb",  "Lu",  "Hf",  "Ta",  "W",   
"Re",  "Os",  "Ir",  "Pt",  "Au",   
"Hg",  "Tl",  "Pb",  "Bi",  "Po",   
"At",  "Rn",  "Fr",  "Ra",  "Ac",   
"Th",  "Pa",  "U",  "Np",  "Pu",   
"Am",  "Cm",  "Bk",  "Cf",  "Es",   
"Fm",  "Md",  "No",  "Lr",  "Rf",   
"Db",  "Sg",  "Bh",  "Hs",  "Mt",   
"Ds",  "Rg",  "Uub",  "Uut",  "Uuq",   
		   "Uup",  "Uuh",  "Uus",  "Uuo"};  
char NUCLEARSYMBOL[119][4]={"NN",  "H",  "HE",  "LI",  "BE",   
"B",  "C",  "N",  "O",  "F",   
"NE",  "NA",  "MG",  "AL",  "SI",   
"P",  "S",  "CL",  "AR",  "K",   
"CA",  "SC",  "TI",  "V",  "CR",   
"MN",  "FE",  "CO",  "NI",  "CU",   
"ZN",  "GA",  "GE",  "AS",  "SE",   
"BR",  "KR",  "RB",  "SR",  "Y",   
"ZR",  "NB",  "MO",  "TC",  "RU",   
"RH",  "PD",  "AG",  "CD",  "IN",   
"SN",  "SB",  "TE",  "I",  "XE",   
"CS",  "BA",  "LA",  "CE",  "PR",   
"ND",  "PM",  "SM",  "EU",  "GD",   
"TB",  "DY",  "HO",  "ER",  "TM",   
"YB",  "LU",  "HF",  "TA",  "W",   
"RE",  "OS",  "IR",  "PT",  "AU",   
"HG",  "TL",  "PB",  "BI",  "PO",   
"AT",  "RN",  "FR",  "RA",  "AC",   
"TH",  "PA",  "U",  "NP",  "PU",   
"AM",  "CM",  "BK",  "CF",  "ES",   
"FM",  "MD",  "NO",  "LR",  "RF",   
"DB",  "SG",  "BH",  "HS",  "MT",   
"DS",  "RG",  "UUB",  "UUT",  "UUQ",   
		   "UUP",  "UUH",  "UUS",  "UUO"}; 
char nuclearsymbol[119][4]={"nn",  "h",  "he",  "li",  "be",   
"b",  "c",  "n",  "o",  "f",   
"ne",  "na",  "mg",  "al",  "si",   
"p",  "s",  "cl",  "ar",  "k",   
"ca",  "sc",  "ti",  "v",  "cr",   
"mn",  "fe",  "co",  "ni",  "cu",   
"zn",  "ga",  "ge",  "as",  "se",   
"br",  "kr",  "rb",  "sr",  "y",   
"zr",  "nb",  "mo",  "tc",  "ru",   
"rh",  "pd",  "ag",  "cd",  "in",   
"sn",  "sb",  "te",  "i",  "xe",   
"cs",  "ba",  "la",  "ce",  "pr",   
"nd",  "pm",  "sm",  "eu",  "gd",   
"tb",  "dy",  "ho",  "er",  "tm",   
"yb",  "lu",  "hf",  "ta",  "w",   
"re",  "os",  "ir",  "pt",  "au",   
"hg",  "tl",  "pb",  "bi",  "po",   
"at",  "rn",  "fr",  "ra",  "ac",   
"th",  "pa",  "u",  "np",  "pu",   
"am",  "cm",  "bk",  "cf",  "es",   
"fm",  "md",  "no",  "lr",  "rf",   
"db",  "sg",  "bh",  "hs",  "mt",   
"ds",  "rg",  "uub",  "uut",  "uuq",   
		   "uup",  "uuh",  "uus",  "uuo"}; 		   
char NuclearName[119][16]= {"neutron",  "hydrogen",  "helium",  "lithium",  "beryllium", 
 "boron",  "carbon",  "nitrogen",  "oxygen",  "fluorine", 
 "neon",  "sodium",  "magnesium",  "aluminum",  "silicon", 
 "phosphorus",  "sulfur",  "chlorine",  "argon",  "potassium", 
 "calcium",  "scandium",  "titanium",  "vanadium",  "chromium", 
 "manganese",  "iron",  "cobalt",  "nickel",  "copper", 
 "zinc",  "gallium",  "germanium",  "arsenic",  "selenium", 
 "bromine",  "krypton",  "rubidium",  "strontium",  "yttrium", 
 "zirconium",  "niobium",  "molybdenum",  "technetium",  "ruthenium", 
 "rhodium",  "palladium",  "silver",  "cadmium",  "indium", 
 "tin",  "antimony",  "tellurium",  "iodine",  "xenon", 
 "cesium",  "barium",  "lanthanum",  "cerium",  "praseodymium", 
 "neodymium",  "promethium",  "samarium",  "europium",  "gadolinium", 
 "terbium",  "dysprosium",  "holmium",  "erbium",  "thulium", 
 "ytterbium",  "lutetium",  "hafnium",  "tantalum",  "tungsten", 
 "rhenium",  "osmium",  "iridium",  "platinum",  "gold", 
 "mercury",  "thallium",  "lead",  "bismuth",  "polonium", 
 "astatine",  "radon",  "francium",  "radium",  "actinium", 
 "thorium",  "protactinium",  "uranium",  "neptunium",  "plutonium", 
 "americium",  "curium",  "berkelium",  "californium",  "einsteinium", 
 "fermium",  "mendelevium",  "nobelium",  "lawrencium",  "rutherfordium", 
 "dubnium",  "seaborgium",  "bohrium",  "hassium",  "meitnerium", 
 "darmstadtium",  "roentgenium",  "ununbium",  "ununtrium",  "ununquadium", 
			    "ununpentium",  "ununhexium",  "ununseptium",  "ununoctium"}; 
/*! \brief Find charge of a nucleus, number before label is allowed, first space is replaced by null-terminated charecter
*/
int NuclearCharge(char *nuclearname) {
char cset[] = " 1234567890";
int sshift = std::strspn (nuclearname,cset); //find length of initial substring
//replace first space after symbol with null character
for (int i=sshift;i<strlen(nuclearname);i++) if (nuclearname[i]==' ') nuclearname[i]='\0';
for (int i=0;i<119;i++) {
//we use length to protect from extra white spaces or other stuff
//int slength=strlen(NuclearSymbol[i]);
if (std::strcmp(nuclearname+sshift,NuclearSymbol[i])==0) return i;
if (std::strcmp(nuclearname+sshift,NUCLEARSYMBOL[i])==0) return i;
if (std::strcmp(nuclearname+sshift,nuclearsymbol[i])==0) return i;
}
//FatalError("Nucleus \""<<nuclearname<<"\" is not known"); 
CriticalError("Nucleus \""<<nuclearname<<"\" is not known");
}
int NuclearNZ(int &N, int &Z, char *nuclearname) {
int A=std::atoi(nuclearname);
Z=NuclearCharge(nuclearname); //shift pointer
N=A-Z;
return A;
}
}
#endif

