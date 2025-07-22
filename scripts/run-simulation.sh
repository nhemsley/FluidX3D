#!/bin/bash

# FluidX3D Simulation Runner Script
# Runs simulation and exports results to a uniquely named directory

set -e  # Exit on error

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNS_DIR="${SCRIPT_DIR}/runs"
BIN_DIR="${SCRIPT_DIR}/bin"
EXPORT_DIR="${BIN_DIR}/export/surface"
HUMAN_HASH_CLI="${SCRIPT_DIR}/human-hash-cli/target/release/human-hash-cli"
SHOREVIEW_DIR="${HOME}/projects/fluidsim/shoreview/test_sequences"

# Default parameters
BATHYMETRY_FILE="assets/sphere.stl"
EXPORT_INTERVAL=1
BUILD_TARGET="Linux-X11"
COPY_TO_SHOREVIEW=true

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --bathymetry)
            BATHYMETRY_FILE="$2"
            shift 2
            ;;
        --export-interval)
            EXPORT_INTERVAL="$2"
            shift 2
            ;;
        --build-target)
            BUILD_TARGET="$2"
            shift 2
            ;;
        --no-shoreview)
            COPY_TO_SHOREVIEW=false
            shift
            ;;
        --help)
            echo "Usage: $0 [options]"
            echo "Options:"
            echo "  --bathymetry FILE         Bathymetry STL file (default: assets/sphere.stl)"
            echo "  --export-interval N       Surface export interval (default: 1)"
            echo "  --build-target TARGET     Build target (default: Linux-X11)"
            echo "  --no-shoreview           Don't copy results to shoreview directory"
            echo "  --help                   Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Ensure human-hash-cli is built
if [ ! -f "$HUMAN_HASH_CLI" ]; then
    echo "Building human-hash-cli..."
    (cd "${SCRIPT_DIR}/human-hash-cli" && cargo build --release)
fi

# Generate unique run identifier
RUN_ID=$($HUMAN_HASH_CLI)
RUN_DIR="${RUNS_DIR}/${RUN_ID}"

echo "=== FluidX3D Simulation Runner ==="
echo "Run ID: ${RUN_ID}"
echo "Output directory: ${RUN_DIR}"
echo "Bathymetry: ${BATHYMETRY_FILE}"
echo "Export interval: ${EXPORT_INTERVAL}"
echo "Build target: ${BUILD_TARGET}"
echo "================================="
echo

# Create runs directory if it doesn't exist
mkdir -p "${RUNS_DIR}"

# Clean up previous exports
echo "Cleaning up previous exports..."
rm -f "${EXPORT_DIR}"/surface_00*
rm -f "${BIN_DIR}/FluidX3D"

# Build FluidX3D
echo "Building FluidX3D..."
make "${BUILD_TARGET}"

# Run simulation
echo "Running simulation..."
"${BIN_DIR}/FluidX3D" \
    --load-bathymetry "${BATHYMETRY_FILE}" \
    --autostart \
    --export-surface-interval="${EXPORT_INTERVAL}"

# Check if surface export was successful
if [ ! -d "${EXPORT_DIR}" ] || [ -z "$(ls -A ${EXPORT_DIR})" ]; then
    echo "Error: No surface files were exported"
    exit 1
fi

# Create run directory and copy results
echo "Copying results to ${RUN_DIR}..."
mkdir -p "${RUN_DIR}"
cp -r "${EXPORT_DIR}" "${RUN_DIR}/"

# Save simulation metadata
cat > "${RUN_DIR}/metadata.json" << EOF
{
    "run_id": "${RUN_ID}",
    "timestamp": "$(date -u +"%Y-%m-%dT%H:%M:%SZ")",
    "bathymetry_file": "${BATHYMETRY_FILE}",
    "export_interval": ${EXPORT_INTERVAL},
    "build_target": "${BUILD_TARGET}",
    "hostname": "$(hostname)",
    "git_commit": "$(git rev-parse HEAD 2>/dev/null || echo 'unknown')",
    "git_branch": "$(git branch --show-current 2>/dev/null || echo 'unknown')"
}
EOF

# Copy to shoreview if requested
if [ "$COPY_TO_SHOREVIEW" = true ] && [ -d "${SHOREVIEW_DIR}" ]; then
    echo "Copying to shoreview directory..."
    rm -rf "${SHOREVIEW_DIR}/surface/"
    cp -r "${EXPORT_DIR}" "${SHOREVIEW_DIR}/"
    echo "Results copied to ${SHOREVIEW_DIR}/surface/"
fi

echo
echo "=== Simulation Complete ==="
echo "Results saved to: ${RUN_DIR}"
echo "Run ID: ${RUN_ID}"

# List exported files
echo
echo "Exported files:"
ls -la "${RUN_DIR}/surface/" | head -10
if [ $(ls "${RUN_DIR}/surface/" | wc -l) -gt 10 ]; then
    echo "... and $(( $(ls "${RUN_DIR}/surface/" | wc -l) - 10 )) more files"
fi

echo
echo "To view in shoreview:"
echo "  cd ~/projects/fluidsim/shoreview"
echo "  # Launch shoreview with the copied surface directory"
