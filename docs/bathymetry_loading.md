# Bathymetry Loading in FluidX3D

This document describes how to use the `--load-bathymetry` command line argument to load seafloor/terrain geometry from STL files for fluid simulations.

## Overview

The bathymetry loading feature allows you to import complex seafloor or terrain geometries from STL files. The simulation domain size is automatically calculated based on the STL file's bounding box, ensuring your geometry fits properly with appropriate padding.

## Usage

```bash
./FluidX3D --load-bathymetry path/to/bathymetry.stl
```

## How It Works

1. **STL Loading**: The STL file is loaded without auto-scaling to preserve the original dimensions
2. **Bounding Box Calculation**: The minimum and maximum coordinates of the geometry are determined
3. **Domain Sizing**: The simulation domain is sized to be 20% larger than the geometry (10% padding on each side)
4. **Resolution Calculation**: Grid resolution is automatically calculated to target ~2GB VRAM usage while maintaining the aspect ratio
5. **Geometry Voxelization**: The STL mesh is voxelized as solid boundary (TYPE_S) in the simulation

## Features

- **Automatic Domain Sizing**: No need to manually specify domain dimensions
- **Aspect Ratio Preservation**: The domain maintains the same aspect ratio as your bathymetry
- **Memory Optimization**: Resolution is automatically calculated to efficiently use available VRAM
- **Wave Generation**: Includes equilibrium boundary conditions for wave generation at y=0

## Example Workflow

1. Create or obtain an STL file of your seafloor/terrain geometry
2. Run the simulation with bathymetry loading:
   ```bash
   ./FluidX3D --load-bathymetry seafloor.stl
   ```
3. The simulation will:
   - Load the STL file
   - Print the bounding box information
   - Calculate appropriate domain size and resolution
   - Set up the fluid simulation with your bathymetry

## Combined with Surface Export

You can combine bathymetry loading with surface export to visualize wave interactions:

```bash
./FluidX3D --load-bathymetry seafloor.stl --export-surface-to ./output/ --export-surface-interval 100
```

## Requirements

- Binary STL file format (standard STL export from most CAD software)
- Required extensions in defines.hpp:
  - `FP16S`
  - `VOLUME_FORCE`
  - `EQUILIBRIUM_BOUNDARIES`
  - `SURFACE` (for free surface simulations)
  - `SURFACE_EXPORT` (if exporting surface data)

## Technical Details

- **Padding**: 10% padding is added on all sides of the bathymetry
- **Target VRAM**: Default target is 2GB, adjustable in the code
- **Boundary Conditions**: 
  - Solid walls on all domain edges
  - Equilibrium boundary at y=0 for wave generation
  - Bathymetry geometry voxelized as solid TYPE_S

## Limitations

- STL file must be in binary format
- Large geometries may require significant VRAM
- The automatic resolution calculation assumes FP16S precision

## Tips

- Ensure your STL file uses consistent units (the simulation assumes dimensionless LBM units)
- For best results, orient your bathymetry so waves propagate in the +y direction
- The water level is set at z = Nz/2 by default