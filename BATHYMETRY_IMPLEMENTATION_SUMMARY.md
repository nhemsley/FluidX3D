# Bathymetry Loading Implementation Summary

## Overview

Implemented a new `--load-bathymetry` command line argument for FluidX3D that allows users to load seafloor or terrain geometry from STL files. The simulation domain is automatically sized based on the STL bounding box, making it easy to simulate waves interacting with complex bathymetry.

## Changes Made

### 1. Modified `src/setup.cpp`

Added bathymetry loading functionality to the breaking waves simulation:

- **Command Line Parsing**: Added code to parse `--load-bathymetry` argument
- **STL Loading**: Loads STL file without auto-scaling to preserve original dimensions
- **Domain Sizing**: Automatically calculates domain size based on STL bounding box with 10% padding
- **Resolution Calculation**: Computes grid resolution to target ~2GB VRAM usage while maintaining aspect ratio
- **Geometry Voxelization**: Voxelizes loaded bathymetry as solid boundaries (TYPE_S)

Key implementation details:
```cpp
// Check for --load-bathymetry command line argument
string bathymetry_file = "";
for(size_t i = 0; i < main_arguments.size(); i++) {
    if(main_arguments[i] == "--load-bathymetry" && i + 1 < main_arguments.size()) {
        bathymetry_file = main_arguments[i + 1];
        break;
    }
}

// Load and process bathymetry if file provided
if(!bathymetry_file.empty()) {
    Mesh* bathymetry_mesh = read_stl(bathymetry_file, 1.0f);
    // Calculate domain size based on bounding box
    // Voxelize as solid geometry
}
```

### 2. Modified `src/setup.hpp`

Added include for `graphics.hpp` to access `main_arguments` global variable for command line parsing.

### 3. Created Documentation

#### `docs/bathymetry_loading.md`
Comprehensive documentation covering:
- Usage instructions
- How the feature works
- Command line examples
- Technical details
- Requirements and limitations
- Tips for best results

#### `examples/bathymetry_example.md`
Practical example guide including:
- Step-by-step usage
- Example workflows
- Creating test bathymetry
- Troubleshooting tips
- Advanced usage options

#### Updated `DOCUMENTATION.md`
Added section about bathymetry loading in the STL loading documentation.

### 4. Created Test Script

`test_bathymetry.sh` - Bash script to test the bathymetry loading functionality with proper error checking and usage examples.

## Usage

Basic usage:
```bash
./FluidX3D --load-bathymetry path/to/bathymetry.stl
```

With surface export:
```bash
./FluidX3D --load-bathymetry seafloor.stl --export-surface-to ./output/ --export-surface-interval 100
```

## Features

1. **Automatic Domain Sizing**: No need to manually calculate domain dimensions
2. **Aspect Ratio Preservation**: Maintains the natural proportions of your bathymetry
3. **Memory Optimization**: Automatically targets efficient VRAM usage
4. **Wave Generation**: Includes equilibrium boundaries for realistic wave generation
5. **Compatibility**: Works with existing surface export functionality

## Technical Details

- **Padding**: 10% padding added on all sides of the bathymetry
- **Target VRAM**: Default 2GB (configurable in code)
- **Water Level**: Set at z = Nz/2 (half domain height)
- **Wave Boundary**: Equilibrium boundary at y=0 for wave generation
- **Boundary Conditions**: Solid walls on all domain edges

## Requirements

Required extensions in `defines.hpp`:
- `FP16S` - Half precision storage
- `VOLUME_FORCE` - For gravity
- `EQUILIBRIUM_BOUNDARIES` - For wave generation
- `SURFACE` - For free surface tracking
- `SURFACE_EXPORT` (optional) - For exporting surface data

## Benefits

1. **Ease of Use**: Single command line argument to load complex geometry
2. **Flexibility**: Works with any binary STL file
3. **Automation**: No manual domain sizing calculations needed
4. **Integration**: Seamlessly integrates with existing FluidX3D features
5. **Visualization**: Compatible with surface export for post-processing

## Future Enhancements

Potential improvements could include:
- Support for multiple STL files (composite bathymetry)
- Automatic water level detection based on geometry
- Variable resolution zones for detailed features
- ASCII STL support
- Terrain texture mapping for visualization