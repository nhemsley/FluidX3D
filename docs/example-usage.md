# FluidX3D Simulation Runner - Example Usage

This document provides examples of using the `run-simulation.sh` script to automate FluidX3D simulations with human-readable run identifiers.

## Prerequisites

1. Build the human-hash-cli tool:
```bash
cd human-hash-cli
cargo build --release
cd ..
```

2. Ensure FluidX3D is properly configured and can build successfully.

## Basic Usage

### Run with Default Settings

```bash
./run-simulation.sh
```

This will:
- Generate a unique run ID (e.g., `autumn-waterfall-flows`)
- Build FluidX3D for Linux-X11
- Run simulation with `assets/sphere.stl` bathymetry
- Export surface every frame
- Save results to `runs/autumn-waterfall-flows/`
- Copy results to shoreview directory

### Run with Custom Bathymetry

```bash
./run-simulation.sh --bathymetry assets/my-terrain.stl
```

### Run with Different Export Interval

```bash
# Export surface every 10 frames
./run-simulation.sh --export-interval 10
```

### Run without Copying to Shoreview

```bash
./run-simulation.sh --no-shoreview
```

### Run with Different Build Target

```bash
# For macOS
./run-simulation.sh --build-target macOS

# For Windows
./run-simulation.sh --build-target Windows
```

## Advanced Examples

### Batch Runs with Different Parameters

```bash
#!/bin/bash
# batch-runs.sh

for interval in 1 5 10 20; do
    echo "Running with export interval $interval"
    ./run-simulation.sh --export-interval $interval
done
```

### Run with Custom Naming

```bash
#!/bin/bash
# Generate a run ID with timestamp
RUN_ID="$(date +%Y%m%d-%H%M%S)-$(./human-hash-cli/target/release/human-hash-cli --words 2)"
mkdir -p "runs/$RUN_ID"

# Run simulation and manually manage output
./bin/FluidX3D --load-bathymetry assets/sphere.stl --autostart --export-surface-interval=1
cp -r bin/export/surface/ "runs/$RUN_ID/"
```

### Compare Multiple Runs

```bash
#!/bin/bash
# compare-runs.sh

# Run simulation with different bathymetries
for bathymetry in sphere.stl cube.stl terrain.stl; do
    ./run-simulation.sh --bathymetry "assets/$bathymetry"
done

# List all runs
echo "Completed runs:"
ls -la runs/
```

## Directory Structure

After running simulations, your directory structure will look like:

```
FluidX3D/
├── runs/
│   ├── autumn-waterfall-flows/
│   │   ├── surface/
│   │   │   ├── surface_000000.stl
│   │   │   ├── surface_000001.stl
│   │   │   └── ...
│   │   └── metadata.json
│   ├── misty-forest-grows/
│   │   ├── surface/
│   │   └── metadata.json
│   └── silent-moon-glows/
│       ├── surface/
│       └── metadata.json
└── run-simulation.sh
```

## Metadata File

Each run includes a `metadata.json` file with information about the simulation:

```json
{
    "run_id": "autumn-waterfall-flows",
    "timestamp": "2024-01-15T10:30:45Z",
    "bathymetry_file": "assets/sphere.stl",
    "export_interval": 1,
    "build_target": "Linux-X11",
    "hostname": "workstation",
    "git_commit": "a1b2c3d4e5f6",
    "git_branch": "main"
}
```

## Viewing Results

### In Shoreview

If shoreview is installed at `~/projects/fluidsim/shoreview/`:

```bash
cd ~/projects/fluidsim/shoreview
# Launch shoreview to view the surface sequence
```

### Direct Access

Access specific run results:

```bash
# List all runs
ls runs/

# View specific run
ls runs/autumn-waterfall-flows/surface/

# Open in mesh viewer
meshlab runs/autumn-waterfall-flows/surface/surface_000000.stl
```

## Troubleshooting

### human-hash-cli not found

```bash
cd human-hash-cli
cargo build --release
```

### No surface files exported

Check that:
1. FluidX3D is built with SURFACE extension enabled
2. The simulation runs long enough to export frames
3. The export directory has write permissions

### Build fails

Ensure all dependencies are installed:
```bash
# For Linux
sudo apt-get install make g++ opencl-headers

# Check available build targets
make help
```

## Tips

1. **Memorable Run IDs**: The human-readable IDs make it easy to remember and discuss specific simulation runs.

2. **Organize by Date**: Create dated subdirectories:
   ```bash
   RUN_DIR="runs/$(date +%Y-%m-%d)/$(./human-hash-cli/target/release/human-hash-cli)"
   ```

3. **Add Notes**: Add a notes.txt file to document specific run parameters:
   ```bash
   echo "Testing new boundary conditions" > "runs/$RUN_ID/notes.txt"
   ```

4. **Clean Old Runs**: Periodically clean up old runs:
   ```bash
   # Remove runs older than 7 days
   find runs/ -type d -mtime +7 -exec rm -rf {} +
   ```
