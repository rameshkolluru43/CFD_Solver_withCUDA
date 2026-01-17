#!/bin/bash
# Google Colab Setup Script for CFD_Solver_withCUDA
# Run this in a Colab notebook cell with: !bash setup_colab.sh

echo "=========================================="
echo "CFD Solver with CUDA - Colab Setup"
echo "=========================================="

# 1. Check CUDA installation
echo -e "\n[1/6] Checking CUDA installation..."
nvcc --version
if [ $? -ne 0 ]; then
    echo "ERROR: CUDA not found!"
    exit 1
fi

# Display GPU info
nvidia-smi

# 2. Update package lists
echo -e "\n[2/6] Updating package lists..."
apt-get update -qq

# 3. Install build tools
echo -e "\n[3/6] Installing build tools..."
apt-get install -y cmake build-essential

# 4. Install JsonCpp
echo -e "\n[4/6] Installing JsonCpp..."
apt-get install -y libjsoncpp-dev

# 5. Install VTK
echo -e "\n[5/6] Installing VTK (this may take a few minutes)..."
apt-get install -y libvtk9-dev vtk9

# Verify installations
echo -e "\n[6/6] Verifying installations..."
cmake --version
g++ --version
nvcc --version

echo -e "\n=========================================="
echo "Dependencies installed successfully!"
echo "=========================================="
echo -e "\nNext steps:"
echo "1. Use CMakeLists_Colab.txt for building"
echo "2. Run: mkdir -p build && cd build"
echo "3. Run: cmake .. -DCMAKE_BUILD_TYPE=Release"
echo "4. Run: make -j$(nproc)"
echo "=========================================="
