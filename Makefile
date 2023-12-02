#define the compilers
CXX=g++
#directory for include
INCDIR=$(COSMO)/include
INCDIR1=$(COSMO)/eigen

#directory to save output programs
BINDIR=$(COSMO)/bin

#this goes for all
CFLAGS= -std=c++17 -fopenmp
INCFLAG=  -I $(INCDIR) -I $(INCDIR1) 
OUTPUTFILE = $(BINDIR)/$@   
.cpp:
	$(CXX) $(CFLAGS) $(INCFLAG)  $@.cpp \
	$(LDFLAGS)  -o $(OUTPUTFILE) 
