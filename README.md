# CFD Solver with CUDA Kernels

A high-performance Computational Fluid Dynamics (CFD) solver featuring GPU acceleration through CUDA kernels. The solver supports both Euler and Navier-Stokes equations with various numerical schemes and boundary conditions for compressible flow simulations.

## 🚀 Features

### Numerical Methods
- **Flux Schemes**: AUSM, Roe, Van Leer, LLF (Local Lax-Friedrichs)
- **Time Integration**: Explicit Runge-Kutta (RK4), TVD-RK3
- **Spatial Discretization**: Second-order accurate with MUSCL reconstruction
- **Limiters**: Van Leer, Minmod, and other slope limiters
- **Gradient Calculation**: Green-Gauss and least squares methods
- **WENO Schemes**: High-order Weighted Essentially Non-Oscillatory methods

### GPU Acceleration
- **CUDA Kernels**: Optimized kernels for flux calculations, gradient computation, time integration
- **Memory Management**: Efficient host-device memory transfers
- **Parallel Execution**: Multi-GPU support with CUDA architectures 6.0-9.0
- **Iterative Solvers**: GPU-accelerated linear algebra operations

### Solver Capabilities
- **Compressible Flow**: Euler and Navier-Stokes equations
- **Boundary Conditions**: Far-field, wall, inlet, outlet conditions
- **Test Cases**: Comprehensive suite of validation cases
- **Grid Support**: Unstructured grids via GMSH format
- **Output Formats**: VTK for visualization

## 📁 Project Structure

```
├── src/                           # Main source files
│   ├── Main.cpp                   # CPU main entry point
│   ├── Main_CUDA.cu              # GPU main entry point
│   ├── Solver.cpp                # Main solver routines
│   ├── Numerical_Method.cpp      # Numerical schemes
│   ├── Boundary_Conditions.cpp   # BC implementations
│   └── ...                       # Additional solver components
├── CUDA_KERNELS/                  # GPU kernel implementations
│   ├── Flux_Calculations_Cuda_Kernels.cu
│   ├── Time_Integration_Cuda_Kernels.cu
│   ├── Gradient_Calculation_Cuda_Kernels.cu
│   ├── Viscous_Flux_Cuda_Kernels.cu
│   ├── Iterative_Solver_Cuda_Kernels.cu
│   └── ...                       # Additional GPU kernels
├── include/                       # Header files
│   ├── definitions.h             # Global definitions
│   ├── Globals.h                 # Global variables
│   ├── Solver.h                  # Solver declarations
│   └── ...                       # Additional headers
├── Test_Cases/                    # Validation test cases
│   ├── Half_Cylinder_Test_Case.cpp
│   ├── Shock_Tube_2D.cpp
│   ├── Flow_Over_Bump.cpp
│   └── ...                       # Additional test cases
├── json_Files/                    # Configuration files
│   ├── Solver_Config.json        # Main solver configuration
│   ├── Boundary_Conditions.json  # BC settings
│   └── ...                       # Test case configurations
├── Grid_Files/                    # Mesh files
├── Gmsh_Grids/                   # GMSH format grids
├── build/                        # Build directory
└── docs/                         # Documentation
```

## 🛠️ Dependencies

### Required
- **CMake** (≥3.16): Build system
- **CUDA Toolkit** (≥11.0): GPU computation
- **C++ Compiler**: C++17 compatible (GCC/Clang)
- **JsonCpp**: JSON configuration parsing
- **VTK** (≥9.4): Visualization output

### Optional
- **Doxygen**: Documentation generation
- **GMSH**: Mesh generation
- **ParaView**: Visualization

## 🔧 Installation

### Prerequisites (macOS with Homebrew)
```bash
# Install dependencies
brew install cmake jsoncpp vtk doxygen

# Install CUDA Toolkit from NVIDIA
# Download from: https://developer.nvidia.com/cuda-downloads
```

### Build Instructions
```bash
# Clone the repository
git clone https://github.com/rameshkolluru43/CFD_Solver_withCUDA.git
cd CFD_Solver_withCUDA

# Create build directory
mkdir -p build
cd build

# Configure with CMake
cmake ..

# Build the project
make -j$(nproc)

# Two executables will be generated:
# - CFD_solver     (CPU version)
# - CFD_solver_gpu (GPU version)
```

## 🚀 Usage

### Running Simulations

#### CPU Version
```bash
./CFD_solver ../json_Files/Solver_Config.json
```

#### GPU Version
```bash
./CFD_solver_gpu ../json_Files/Solver_Config.json
```

### Configuration

The solver uses JSON configuration files to define simulation parameters:

#### Main Configuration (`Solver_Config.json`)
```json
{
    "TestCase": {
        "Test_Case": 1,
        "Test_Case_Name": "Flow_Over_Cylinder",
        "Test_Case_Json": "../json_Files/Half_Cylinder.json"
    },
    "Simulation": {
        "Total_Iterations": 100000,
        "CFL": 0.1,
        "Is_Time_Dependent": false
    },
    "Solver": {
        "Solver_Type": 0,
        "Is_Viscous": false,
        "Flux_Type": 2,
        "Is_Second_Order": true
    }
}
```

### Available Test Cases

1. **Half Cylinder**: Flow over a half cylinder
2. **Shock Tube**: 1D shock tube problem
3. **Flow Over Bump**: Transonic flow over a bump
4. **Forward Facing Step**: Supersonic flow over a step
5. **Shock Reflection**: 2D shock reflection
6. **SWBLI**: Shock-Wave Boundary Layer Interaction
7. **Scramjet Inlet**: Hypersonic inlet flow
8. **And more...**

## 🧪 Validation & Testing

The solver includes comprehensive test cases and validation framework:

### Grid Processing Tests ✅
- **VTK Grid Reading**: Comprehensive testing with `test_read_vtk_grid_integrated.cpp`
- **Test Grid**: `Ramp_15o_52_18.vtk` (867 points, 800 quadrilateral cells)
- **Geometric Validation**: Face normals, cell areas, centroids with real CFD solver functions
- **Boundary Data Structures**: Ghost cell construction, boundary classification
- **Data Structure Integrity**: Co-Volume cells, cell-ghost connectivity validation

### Test Results Summary
- **Physical Cells**: 800 with proper geometric properties
- **Ghost Cells**: 131 for boundary conditions  
- **Boundary Classification**: 16 inlet, 16 exit, 99 wall boundaries
- **Geometric Quality**: Area ratio 1.86 (well-conditioned grid)
- **Face Normals**: Unit vectors with proper outward orientation
- **Area Conservation**: Validated with `Check_Cells()` function

### Flow Validation Cases
- Shock tube problems for accuracy verification
- Boundary layer flows for viscous validation
- Complex geometries for robustness testing
- Comparison with analytical solutions where available

For detailed testing information, see [TEST_SUMMARY.md](TEST_SUMMARY.md).

## 📊 GPU Performance

The CUDA implementation provides significant speedup over CPU execution:
- **Memory Bandwidth**: Optimized memory access patterns
- **Kernel Optimization**: Coalesced memory access, shared memory usage
- **Multi-GPU Support**: Scalable to multiple GPUs
- **Architecture Support**: CUDA compute capability 6.0-9.0

### Supported CUDA Architectures
- Pascal (6.0, 6.1)
- Volta (7.0)
- Turing (7.5)
- Ampere (8.0, 8.6)
- Ada Lovelace (8.9)
- Hopper (9.0)

## 📖 Documentation

### Generate Documentation
```bash
# Using Doxygen
doxygen Doxyfile_Cleaned

# Documentation will be generated in:
# - docs/html/     (HTML format)
# - docs/latex/    (LaTeX format)
```

### Key Classes and Functions
- `Solver`: Main solver class
- `Cell`: Computational cell representation
- `Face`: Face/interface handling
- `Flux`: Numerical flux calculations
- `Boundary_Conditions`: Boundary condition implementations

## 🔬 Development

### Code Structure
- **Object-Oriented Design**: Modular, maintainable code
- **CUDA Integration**: Seamless CPU-GPU data transfer
- **Memory Management**: Efficient memory usage patterns
- **Error Handling**: Comprehensive error checking

### Adding New Features
1. **New Numerical Schemes**: Add to `src/Numerical_Method.cpp`
2. **New Boundary Conditions**: Extend `src/Boundary_Conditions.cpp`
3. **New CUDA Kernels**: Add to appropriate `CUDA_KERNELS/*.cu` file
4. **New Test Cases**: Create in `Test_Cases/` directory

## 📄 License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## 👥 Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## 📞 Contact

**Ramesh Kolluru** - rameshkolluru43@gmail.com

Project Link: [https://github.com/rameshkolluru43/CFD_Solver_withCUDA](https://github.com/rameshkolluru43/CFD_Solver_withCUDA)

## 🙏 Acknowledgments

- CUDA programming community for optimization techniques
- CFD community for validation test cases
- Open source contributors for tools and libraries used

## 📚 References

- Computational Fluid Dynamics literature
- CUDA programming best practices
- Numerical methods for hyperbolic equations

---

**Note**: This solver is actively developed for research purposes. For production use, thorough validation is recommended for your specific applications.
