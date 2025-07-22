#!/bin/bash

# Script to analyze bathymetry simulation frames

echo "Analyzing bathymetry simulation frames..."
echo "========================================="

# Check if frames directory exists
FRAMES_DIR="../tmp/bathymetry_frames"
if [ ! -d "$FRAMES_DIR" ]; then
    echo "Error: Frames directory not found at $FRAMES_DIR"
    exit 1
fi

# Count frames
FRAME_COUNT=$(ls -1 $FRAMES_DIR/*.png 2>/dev/null | wc -l)
echo "Total frames found: $FRAME_COUNT"

if [ $FRAME_COUNT -eq 0 ]; then
    echo "No PNG frames found in $FRAMES_DIR"
    exit 1
fi

# Get first and last frame
FIRST_FRAME=$(ls -1 $FRAMES_DIR/*.png | head -1)
LAST_FRAME=$(ls -1 $FRAMES_DIR/*.png | tail -1)

echo "First frame: $(basename $FIRST_FRAME)"
echo "Last frame: $(basename $LAST_FRAME)"
echo ""

# Create comparison directory
COMPARE_DIR="../tmp/frame_comparison"
mkdir -p $COMPARE_DIR

# Copy key frames for analysis
echo "Copying key frames for analysis..."
cp "$FIRST_FRAME" "$COMPARE_DIR/01_first_frame.png"
cp "$LAST_FRAME" "$COMPARE_DIR/99_last_frame.png"

# Copy some intermediate frames
for timestep in 0500 1000 1500 2000 2500 3000 3500 4000 4500; do
    FRAME_FILE="$FRAMES_DIR/image-00000${timestep}.png"
    if [ -f "$FRAME_FILE" ]; then
        cp "$FRAME_FILE" "$COMPARE_DIR/frame_${timestep}.png"
        echo "  Copied frame at timestep $timestep"
    fi
done

echo ""
echo "Key frames copied to: $COMPARE_DIR"
echo ""

# Check if ImageMagick is installed for frame analysis
if command -v identify &> /dev/null; then
    echo "Frame properties (using ImageMagick):"
    identify -format "  Resolution: %wx%h\n  Format: %m\n  File size: %b\n" "$FIRST_FRAME"
else
    echo "Install ImageMagick for detailed frame analysis"
fi

echo ""
echo "Frame file sizes:"
ls -lh $COMPARE_DIR/*.png | awk '{print "  " $9 ": " $5}'

echo ""
echo "Done! Check $COMPARE_DIR for frame comparisons"
