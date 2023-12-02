# This script must be sourced reference should not have ".." e.g. ./../path
#--------------DETERMINE LOCATION OF SOURCED SCRIPT------------
MYLOCATION=$0
#the following commad finds the directory of sourced script
#if below protects from filename being used instead of directory
if [ -f "$MYLOCATION" ]
then
 FULLPATH=$(pwd)
else 
FULLPATH=$(echo "$MYLOCATION" | sed  "s|\.|$(pwd)|")
fi

 
#----------LOCATION IS DETERMINED------------------

#----------SET PATH VARIABLE DEPENDING ON ARGUMENT-----------
if [ -d "$1"  ]; 
then 
export COSMO=$1
#use parameter to set path
else 
export COSMO=$FULLPATH
fi

if [ -d "$COSMO"  ]
then 
#ls -d $PWD/*
echo "Shell Model Path:" $COSMO
#make sure we bave bin path
if [ ! -d "$COSMO/bin" ]; then
  mkdir -p "$COSMO/bin";
fi
#export path 
export PATH="$COSMO/bin:$PATH"
#export
if [ ! -d "$COSMO/interactions" ]; then
  echo "--------Problem---------"
  echo "Interactions are expected in $COSMO/interactions"
  echo "but were not found, most defaults will not work"
  echo "SMINTPATH  variable is not set"
  echo "-------------------------"
else
SMINTPATH=$(ls -d $COSMO/interactions/*/ | awk '{printf ":"$0}' | sed '1s/^.//')
export SMINTPATH="$COSMO/interactions:$SMINTPATH"
fi
export JISPPATH="$COSMO/JISPintFiles"
if [ -f "$COSMO/Makefile" ]; then  
export MAKEFILES="$COSMO/Makefile"
echo "Default Makefile:" $MAKEFILES
else
echo "Default Makefile is unchanged"
#advise on makefile
	if [ -f "$COSMO/Makefile_template" ]; then
	echo "--------Notice---------"
	echo "If you are running for the first time consider "
	echo "copying Makefile_template to Makefile "
	echo "and setting CXX to and flags to match your system"
	echo "-----------------------"
	fi
fi

else
echo "No valid direcoty is found"
echo "Either run script by providing a full path to it or use an argument" 
echo "Argument: $1"
fi

