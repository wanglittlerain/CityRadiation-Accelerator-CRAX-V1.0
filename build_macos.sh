#!/bin/bash

# macOS build script for CRAX
# This script builds the project with macOS-specific configurations

set -e  # Exit on any error

echo "Building CRAX for macOS..."

# Create necessary directories
mkdir -p libs bin

# Clean previous build
if [ -d "build" ]; then
    echo "Cleaning previous build..."
    rm -rf build/*
fi

# Configure with CMake
echo "Configuring with CMake..."
cmake -B build -S .

# Build the project
echo "Building..."
cd build && make

echo "Build completed successfully!"
echo "Executables available in bin/:"
ls -l ../bin/

echo ""
echo "To run the programs:"
echo "  ./bin/radiation [args]"
echo "  ./bin/mesh-generation [args]"
