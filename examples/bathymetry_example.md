# Bathymetry Loading Example

This example demonstrates how to use the `--load-bathymetry` feature in FluidX3D to load seafloor or terrain geometry from STL files for breaking waves simulations.

## Prerequisites

1. Build FluidX3D with the following extensions enabled in `defines.hpp`:
   - `FP16S`
   - `VOLUME_FORCE`
   - `EQUILIBRIUM_BOUNDARIES`
   - `SURFACE`
   - `INTERACTIVE_GRAPHICS` (optional, for visualization)
   - `SURFACE_EXPORT` (optional, for exporting surface data)

2. Have a binary STL file of your bathymetry/terrain ready

## Basic Usage

### Command Line

```bash
./bin/FluidX3D --load-bathymetry path/to/your/bathymetry.stl
```

### What Happens

1. The STL file is loaded and its bounding box is calculated
2. The simulation domain is automatically sized to be 20% larger than the bathymetry (10% padding on each side)
3. Grid resolution is calculated to efficiently use ~2GB of VRAM while maintaining the aspect ratio
4. The bathymetry is voxelized as solid boundaries in the simulation
5. Water is initialized at half the domain height with hydrostatic pressure
6. Waves are generated at the y=0 boundary

## Example with Surface Export

To visualize the waves interacting with your bathymetry:

```bash
./bin/FluidX3D --load-bathymetry seafloor.stl --export-surface-to ./output/ --export-surface-interval 100
```

This will:
- Load the bathymetry from `seafloor.stl`
- Export the fluid surface as STL files every 100 timesteps
- Save the surface files to the `./output/` directory

## Creating Test Bathymetry

You can create simple bathymetry STL files using any 3D modeling software. Here's a basic workflow:

1. **Flat Seafloor with Ridge**:
   - Create a large flat plane
   - Add a ridge or mountain in the middle
   - Export as binary STL

2. **Sloped Beach**:
   - Create an inclined plane
   - Add some rocks or features
   - Export as binary STL

3. **Complex Terrain**:
   - Use terrain generation tools
   - Import height maps
   - Export as binary STL

## Example Output

When you run the simulation with bathymetry loading, you'll see output like:

```
Loading bathymetry from: seafloor.stl
Bathymetry bounding box:
  Min: -50.0, -100.0, -20.0
  Max: 50.0, 100.0, 5.0
  Size: 100.0 x 200.0 x 25.0
Domain resolution: 256 x 512 x 64
```

## Tips

1. **Units**: The simulation uses dimensionless LBM units. Scale your STL appropriately.

2. **Orientation**: Orient your bathymetry so waves will propagate in the +y direction.

3. **Water Level**: By default, water fills up to half the domain height (z = Nz/2).

4. **Memory Usage**: Large bathymetries will require more VRAM. Adjust the target VRAM in the code if needed.

5. **Wave Generation**: Waves are generated at the y=0 boundary using an equilibrium boundary condition with sinusoidal velocity.

## Troubleshooting

- **"File not found"**: Ensure the STL file path is correct
- **Out of memory**: Reduce the domain size by modifying the target VRAM in the code
- **No waves**: Check that your bathymetry doesn't block the y=0 boundary
- **Simulation crashes**: Ensure all required extensions are enabled in `defines.hpp`

## Advanced Usage

You can modify the wave parameters in the `main_setup()` function:
- `f`: Gravity force (default: 0.001f)
- `u`: Peak velocity of wave generation (default: 0.12f)
- `frequency`: Wave frequency (default: 0.0007f)

## Example STL Files

You can find example STL files for testing at:
- Simple geometries: https://www.thingiverse.com/
- Terrain data: Convert DEM/height maps to STL using online tools
- Ocean bathymetry: Process NOAA bathymetry data into STL format