#!/bin/bash

# Test script for FluidX3D mesh streaming with seaview server
# Runs seaview on a unique port for testing mesh streaming functionality

set -e  # Exit on error

# Configuration
DEFAULT_PORT=8089
SEAVIEW_BIN="seaview"
LOG_FILE="/tmp/seaview_test_streaming.log"

# Parse command line arguments
PORT=$DEFAULT_PORT
while [[ $# -gt 0 ]]; do
    case $1 in
        --port)
            PORT="$2"
            shift 2
            ;;
        --help)
            echo "Usage: $0 [options]"
            echo "Options:"
            echo "  --port PORT    Port to run seaview server on (default: $DEFAULT_PORT)"
            echo "  --help         Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

echo "=== Seaview Test Streaming Server ==="
echo "Port: $PORT"
echo "Log file: $LOG_FILE"
echo "===================================="
echo

# Check if seaview is installed
if ! command -v "$SEAVIEW_BIN" &> /dev/null; then
    echo "Error: seaview not found in PATH"
    echo "Please ensure seaview is installed and available in your PATH"
    exit 1
fi

# Kill any existing seaview processes on this port
echo "Checking for existing seaview processes on port $PORT..."
if lsof -i :$PORT &> /dev/null; then
    echo "Found existing process on port $PORT, killing it..."
    lsof -ti :$PORT | xargs kill -9 2>/dev/null || true
    sleep 1
fi

# Start seaview server
echo "Starting seaview server on port $PORT..."
echo "Server will accept mesh streaming connections from FluidX3D"
echo
echo "To test mesh streaming:"
echo "1. Run FluidX3D with mesh streaming enabled"
echo "2. FluidX3D will attempt to connect to localhost:$PORT"
echo "3. Check the log file for connection attempts: tail -f $LOG_FILE"
echo
echo "Press Ctrl+C to stop the server"
echo

# Run seaview with specified port, redirecting output to log file
$SEAVIEW_BIN --port $PORT 2>&1 | tee "$LOG_FILE"
