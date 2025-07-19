# Surface Export Feature Documentation

## Overview

The surface export feature in FluidX3D allows you to export the fluid surface mesh to STL files during simulation. This is useful for:
- Post-processing in external software (Blender, ParaView, etc.)
- 3D printing fluid simulation results
- Creating animations of fluid surfaces
- Scientific visualization and analysis

## Requirements

To use the surface export feature, you must enable the following extensions in `defines.hpp`:

```cpp
#define SURFACE         // Required: enables free surface LBM
#define EXPORT_SURFACE  // Required: enables surface mesh export from GPU
#define SURFACE_EXPORT  // Required: enables STL file writing functionality
```

## Usage

### Method 1: Automatic Export with `run_with_surface_export()`

The simplest way to export surfaces is to use the `run_with_surface_export()` function instead of `lbm.run()`:

```cpp
void main_setup() {
    // Create and configure your LBM simulation
    LBM lbm(128u, 256u, 256u, 0.005f, 0.0f, 0.0f, -0.0002f, 0.0001f);
    
    // Set up initial conditions (fluid, boundaries, etc.)
    // ...
    
    // Run simulation with automatic surface export
#ifdef SURFACE_EXPORT
    run_with_surface_export(&lbm);  // Runs indefinitely with periodic export
    // or
    run_with_surface_export(&lbm, 5000u);  // Run for 5000 timesteps
#else
    lbm.run();
#endif
}
```

### Method 2: Manual Export in Custom Loop

For more control over when surfaces are exported:

```cpp
void main_setup() {
    LBM lbm(128u, 128u, 128u, 0.01f);
    const uint lbm_T = 10000u;
    
    // Set up simulation...
    
    lbm.run(0u);  // Initialize
    
    while(lbm.get_t() < lbm_T) {
        lbm.run(100u);  // Run 100 steps
        
#ifdef SURFACE_EXPORT
        // Export at specific intervals or conditions
        if(lbm.get_t() % 500u == 0u) {
            export_surface_frame(&lbm, lbm.get_t());
        }
#endif
    }
}
```

## Command Line Arguments

The surface export feature is controlled via command line arguments:

### Required Arguments

- `--export-surface-to <directory>` - Directory where STL files will be saved
  - Example: `--export-surface-to ./output/`
  - The directory will be created if it doesn't exist

### Optional Arguments

- `--export-surface-interval <N>` - Export every N timesteps (default: 100)
  - Example: `--export-surface-interval 50`

- `--export-surface-ascii` - Use ASCII STL format instead of binary (default: binary)
  - ASCII files are larger but human-readable
  - Binary files are more compact and faster to write

## Example Command

```bash
./FluidX3D --export-surface-to ./stl_output/ --export-surface-interval 100

# With ASCII format
./FluidX3D --export-surface-to ./stl_output/ --export-surface-interval 100 --export-surface-ascii
```

## Output Files

Surface meshes are saved as STL files with the naming convention:
```
surface_XXXXXX.stl
```

Where `XXXXXX` is the zero-padded timestep number. For example:
- `surface_000000.stl` - Initial state
- `surface_000100.stl` - Timestep 100
- `surface_000200.stl` - Timestep 200

## Performance Considerations

1. **Export Frequency**: Exporting too frequently can impact performance. Choose an appropriate interval based on your needs.

2. **File Size**: 
   - Binary STL files are typically 50 bytes per triangle plus 84 bytes header
   - ASCII STL files can be 5-10x larger
   - A complex fluid surface might have 100k-1M triangles

3. **Disk Space**: Ensure you have sufficient disk space for the exported files.

4. **GPU Memory**: The surface export requires additional GPU memory to store triangle vertices.

## Troubleshooting

### No Output Files
- Ensure `--export-surface-to` directory is specified
- Check that SURFACE, EXPORT_SURFACE, and SURFACE_EXPORT are all defined
- Verify the directory has write permissions

### Empty STL Files
- Make sure your simulation contains fluid (TYPE_F cells)
- The surface extraction only works where fluid interfaces exist

### Out of Memory
- The surface buffer automatically resizes, but very complex surfaces may exceed GPU memory
- Try reducing the simulation resolution or export frequency

## Example Setups

FluidX3D includes example setups demonstrating surface export:

1. **Dam Break with Surface Export** (in `setup.cpp`)
   ```cpp
   // Uncomment and run with:
   // ./FluidX3D --export-surface-to ./output/ --export-surface-interval 100
   ```

2. **Custom Export Loop** (in `setup.cpp`)
   - Shows manual control over export timing
   - Combines surface export with graphics rendering

## Visualization

The exported STL files can be visualized in:
- **Blender**: Import via File → Import → STL
- **ParaView**: Open directly, apply filters as needed
- **MeshLab**: Good for quick viewing and basic processing
- **CloudCompare**: Useful for point cloud and mesh analysis

## Tips

1. For animations, import the STL sequence into Blender using the "Import Images as Planes" addon or similar tools
2. Use binary format for production runs to save disk space
3. Export less frequently at the beginning if the simulation takes time to develop
4. Consider post-processing to reduce mesh complexity if needed