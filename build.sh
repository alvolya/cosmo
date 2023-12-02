#!/bin/bash

# Step 1: Get the path to the script's location and set COSMO variable
COSMO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Set the C++ compiler to g++
CXX=g++

# Step 2: Check if $COSMO/bin directory exists, if not, create it
if [ ! -d "$COSMO/bin" ]; then
    mkdir "$COSMO/bin"
fi

# Step 3: Check if $COSMO/eigen exists, if not, prompt user to install it
if [ ! -d "$COSMO/eigen" ]; then
    read -p "Eigen library not found. Would you like to install it? [Y/n] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]
    then
        (cd "$COSMO" && git clone https://gitlab.com/libeigen/eigen.git)
    else
        echo "Please install Eigen library in $COSMO/eigen and rerun the script."
        exit 1
    fi
fi

# Step 4: Compile all .cpp files, excluding those in the eigen directory
find "$COSMO" -type f -name '*.cpp' ! -path "$COSMO/eigen/*" ! -path "$COSMO/include/*" | while read -r file; do
    echo "Compiling $file..."
    $CXX -std=c++17 -fopenmp -I"$COSMO/include" -I"$COSMO/eigen" "$file" -o "$COSMO/bin/$(basename "${file%.*}")"
done
