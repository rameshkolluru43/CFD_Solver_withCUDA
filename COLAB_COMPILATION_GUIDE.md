# CFD Solver with CUDA - Google Colab Compilation Guide

This guide explains how to compile and run the CFD Solver with CUDA on Google Colab.

## Prerequisites

Google Colab provides:
- ✅ CUDA Toolkit (usually 11.x or 12.x)
- ✅ NVIDIA GPU (Tesla T4, P100, or V100)
- ✅ Ubuntu Linux environment

## Step 1: Clone Repository

```python
# In a Colab notebook cell
!git clone https://github.com/rameshkolluru43/CFD_Solver_withCUDA.git
%cd CFD_Solver_withCUDA
```

## Step 2: Install Dependencies

### Option A: Use the setup script
```python
!chmod +x setup_colab.sh
!bash setup_colab.sh
```

### Option B: Manual installation
```python
# Check CUDA
!nvcc --version
!nvidia-smi

# Install dependencies
!apt-get update
!apt-get install -y cmake build-essential libjsoncpp-dev libvtk9-dev vtk9
```

## Step 3: Configure CMake

```python
# Copy Colab-specific CMakeLists.txt
!cp CMakeLists_Colab.txt CMakeLists.txt

# OR manually update the original CMakeLists.txt
# (see modifications section below)
```

## Step 4: Build the Project

```python
# Create build directory
!mkdir -p build
%cd build

# Configure with CMake
!cmake .. -DCMAKE_BUILD_TYPE=Release

# Build (use all available cores)
!make -j$(nproc)

# Check if executables were created
!ls -lh CFD_solver_gpu
```

## Step 5: Run Simulations

```python
# Go back to project root
%cd ..

# Run GPU version
!./build/CFD_solver_gpu ./json_Files/Solver_Config.json
```

## Common Issues and Solutions

### Issue 1: CUDA Architecture Mismatch
**Error:** `No kernel image is available for execution on the device`

**Solution:** Check your GPU compute capability:
```python
!nvidia-smi --query-gpu=compute_cap --format=csv,noheader
```

Then modify `CMakeLists_Colab.txt` line with `CUDA_ARCHITECTURES`:
- T4 GPU: `"75"`
- P100 GPU: `"60"`
- V100 GPU: `"70"`
- A100 GPU: `"80"`

### Issue 2: VTK Library Not Found
**Solution:**
```python
# Try installing VTK from pip instead
!pip install vtk

# Or use older VTK version
!apt-get install -y libvtk7-dev
```

### Issue 3: JsonCpp Headers Not Found
**Solution:**
```python
# Install from apt
!apt-get install -y libjsoncpp-dev

# Or verify include path
!find /usr/include -name "json.h" 2>/dev/null | grep jsoncpp
```

### Issue 4: Out of Memory
Colab GPUs have limited memory (15-16GB for T4). If you run out of memory:
```python
# Reduce problem size in json configuration
# Or use CPU version for testing
!./build/CFD_solver ./json_Files/Solver_Config.json
```

## Complete Colab Notebook Example

```python
# Cell 1: Setup and clone
!git clone https://github.com/rameshkolluru43/CFD_Solver_withCUDA.git
%cd CFD_Solver_withCUDA

# Cell 2: Check GPU
!nvidia-smi
!nvcc --version

# Cell 3: Install dependencies
!apt-get update -qq
!apt-get install -y cmake build-essential libjsoncpp-dev libvtk9-dev

# Cell 4: Configure for Colab
!cp CMakeLists_Colab.txt CMakeLists.txt

# Cell 5: Build
!mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && make -j$(nproc)

# Cell 6: Run simulation
!./build/CFD_solver_gpu ./json_Files/Solver_Config.json

# Cell 7: Visualize results (if VTK output is generated)
# Download output files or use in-notebook visualization
from google.colab import files
!zip -r results.zip *.vtk
files.download('results.zip')
```

## Key Differences from Local Compilation

| Aspect | Local (macOS) | Google Colab |
|--------|---------------|--------------|
| OS | macOS | Ubuntu Linux |
| Package Manager | Homebrew | apt-get |
| CUDA Path | `/usr/local/cuda` | `/usr/local/cuda` |
| Include Paths | `/opt/homebrew/include` | `/usr/include` |
| Library Paths | `/opt/homebrew/lib` | `/usr/lib` |
| GPU Architecture | Multiple (60-90) | Specific (75 for T4) |
| VTK Version | 9.4 | 9.x (varies) |

## Optimization Tips for Colab

1. **Persistent Storage**: Mount Google Drive to save build artifacts
   ```python
   from google.colab import drive
   drive.mount('/content/drive')
   ```

2. **GPU Runtime**: Enable GPU in Runtime → Change runtime type → GPU

3. **Build Caching**: Save the build directory to avoid recompiling
   ```python
   !tar -czf build_cache.tar.gz build/
   # Copy to Drive for next session
   ```

4. **Memory Management**: Monitor GPU memory usage
   ```python
   !watch -n 1 nvidia-smi  # Run in background
   ```

## Testing the Installation

Run a quick test to verify everything works:
```python
# Simple compilation test
!cd build && ./CFD_solver_gpu --help

# Quick simulation test with small grid
# Modify json config for smaller problem if needed
```

## Troubleshooting Commands

```python
# Check CMake configuration
!cat build/CMakeCache.txt | grep CUDA

# Check linked libraries
!ldd build/CFD_solver_gpu

# Find library locations
!ldconfig -p | grep -E 'jsoncpp|vtk|cuda'

# Verify CUDA samples can compile
!cd /usr/local/cuda/samples/1_Utilities/deviceQuery && make && ./deviceQuery
```

## Resources

- [Google Colab CUDA Guide](https://colab.research.google.com/notebooks/gpu.ipynb)
- [NVIDIA CUDA Documentation](https://docs.nvidia.com/cuda/)
- [Project GitHub Repository](https://github.com/rameshkolluru43/CFD_Solver_withCUDA)
