#!/bin/bash

# Set the path to your NetCDF source directory
NETCDF_DIR="$1"

# Check if the specified directory exists
if [ ! -d "$NETCDF_DIR" ]; then
    echo "Directory $NETCDF_DIR does not exist."
    exit 1
fi

# Navigate to the NetCDF source directory
cd "$NETCDF_DIR" || exit

# Path to the ncvalues.cpp file
NCVALUES_FILE="cxx/ncvalues.cpp"

# Check if the ncvalues.cpp file exists
if [ ! -f "$NCVALUES_FILE" ]; then
    echo "$NCVALUES_FILE not found."
    exit 1
fi

# Add #include <cstring> if it is not already present
if ! grep -q '#include <cstring>' "$NCVALUES_FILE"; then
    echo "Adding #include <cstring> to $NCVALUES_FILE"
    sed -i 's|#include "ncvalues.h"|#include "ncvalues.h"\n#include <cstring>|' "$NCVALUES_FILE"
else
    echo "#include <cstring> already present in $NCVALUES_FILE"
fi

