#!/bin/bash

# Test script for bathymetry loading functionality in FluidX3D
# This script demonstrates how to use the --load-bathymetry command line argument

echo "============================================"
echo "FluidX3D Bathymetry Loading Test"
echo "============================================"
echo ""

# Check if FluidX3D binary exists
if [ ! -f "bin/FluidX3D" ]; then
    echo "Error: FluidX3D binary not found. Please run ./make.sh first."
    exit 1
fi

# Function to show usage
show_usage() {
    echo "Usage: $0 <path_to_bathymetry.stl>"
    echo ""
    echo "Example:"
    echo "  $0 ../stl/seafloor.stl"
    echo ""
    echo "Additional options can be combined:"
    echo "  $0 ../stl/seafloor.stl --export-surface"
    echo ""
    echo "The domain size will be automatically calculated based on the STL bounding box."
}

# Check if STL file path is provided
if [ $# -lt 1 ]; then
    echo "Error: No STL file specified."
    echo ""
    show_usage
    exit 1
fi

STL_FILE="$1"

# Check if STL file exists
if [ ! -f "$STL_FILE" ]; then
    echo "Error: STL file not found: $STL_FILE"
    exit 1
fi

echo "Loading bathymetry from: $STL_FILE"
echo ""

# Prepare command with bathymetry loading
CMD="./bin/FluidX3D --load-bathymetry \"$STL_FILE\""

# Check if surface export is requested
if [[ "$*" == *"--export-surface"* ]]; then
    # Create export directory if it doesn't exist
    EXPORT_DIR="export/bathymetry_test/"
    mkdir -p "$EXPORT_DIR"

    CMD="$CMD --export-surface-to \"$EXPORT_DIR\" --export-surface-interval 100"
    echo "Surface export enabled to: $EXPORT_DIR"
fi

echo ""
echo "Running command:"
echo "$CMD"
echo ""
echo "============================================"
echo ""

# Run FluidX3D with bathymetry loading
eval $CMD
