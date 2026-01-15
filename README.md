# CFD Solver with CUDA Kernels

A high-performance Computational Fluid Dynamics (CFD) solver featuring GPU acceleration through CUDA kernels. The solver supports both Euler and Navier-Stokes equations with various numerical schemes and boundary conditions for compressible flow simulations.

## 🚀 Features

### Numerical Methods

#### **Advanced Flux Computation Schemes** 🚀
- **Van Leer Flux Vector Splitting**: Complete implementation with Mach number-based splitting for exact contact preservation
- **Roe Approximate Riemann Solver**: 
  - **First-Order**: Enhanced with entropy fix, comprehensive error checking, and boundary condition handling
  - **Second-Order**: TVD implementation with slope limiting and high-resolution shock capturing
- **AUSM (Advection Upstream Splitting Method)**: Robust flux computation for all Mach number regimes
- **LLF (Local Lax-Friedrichs)**: Simple and robust flux approximation

#### **High-Order Methods**
- **Time Integration**: Explicit Runge-Kutta (RK4), TVD-RK3 with optimal stability properties
- **Spatial Discretization**: Second-order accurate with MUSCL reconstruction and limiter integration
- **Slope Limiters**: Van Leer, Minmod, Superbee with TVD properties for shock capturing
- **Gradient Calculation**: Green-Gauss and weighted least squares methods
- **WENO Schemes**: High-order Weighted Essentially Non-Oscillatory methods for complex flows

### GPU Acceleration
- **CUDA Kernels**: Optimized kernels for flux calculations, gradient computation, time integration
- **Memory Management**: Efficient host-device memory transfers
- **Parallel Execution**: Multi-GPU support with CUDA architectures 6.0-9.0
- **Iterative Solvers**: GPU-accelerated linear algebra operations

### Solver Capabilities

#### **Flow Physics** 🌊
- **Compressible Flow**: Full Euler and Navier-Stokes equations with thermodynamic consistency
- **Flow Regimes**: Subsonic, transonic, supersonic, and hypersonic flow support
- **Shock Capturing**: Advanced shock-capturing with entropy-satisfying schemes
- **Contact Preservation**: Exact contact discontinuity preservation with specialized flux methods

#### **Computational Framework** 🔧
- **Boundary Conditions**: Far-field, wall, inlet, outlet, symmetry with robust treatment
- **Grid Support**: Unstructured grids via GMSH format with automatic boundary detection
- **Test Cases**: Comprehensive validation suite including shock tubes, cylinder flows, and complex geometries
- **Output Formats**: VTK for ParaView visualization with field variable export
- **Error Handling**: Production-ready error checking and graceful failure recovery

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
    ├── Van_Leer_Flux_Implementation.md       # Van Leer flux technical documentation
    ├── ROE_2O_Implementation.md               # Second-order Roe scheme documentation  
    ├── Enhanced_ROE_First_Order.md           # Enhanced first-order Roe documentation
    ├── AUSM_Flux_Implementation.md           # AUSM flux scheme documentation
    └── *_Completion_Summary.md                # Implementation completion summaries
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

## � Flux Scheme Enhancements

### **Production-Ready Implementations** ✅

#### **Van Leer Flux Vector Splitting**
- **Complete Implementation**: Mach number-based flux splitting with exact contact preservation
- **Performance**: ~50 floating point operations per face with optimal computational efficiency
- **Flow Regimes**: Excellent performance across subsonic, transonic, and supersonic flows
- **Documentation**: Comprehensive 4000+ word technical guide with mathematical framework

#### **Enhanced Roe Approximate Riemann Solver**
- **First-Order Enhancement**: 
  - Entropy fix for sonic points preventing expansion shocks
  - Comprehensive error checking and state validation
  - Robust boundary condition handling with graceful fallbacks
  - Production-ready reliability with industrial-grade error handling
- **Second-Order Implementation**:
  - TVD slope limiting with multiple limiter options (Van Leer, Minmod, Superbee)
  - High-resolution shock capturing without spurious oscillations
  - Automatic limiter selection based on local flow conditions
  - Mathematical rigor with complete eigenvalue-eigenvector decomposition

#### **AUSM Flux Scheme**
- **All-Speed Capability**: Robust performance from incompressible to hypersonic regimes
- **Mass Flux Splitting**: Advanced upwind splitting for momentum and energy equations
- **Pressure Correction**: Proper treatment of pressure and acoustic waves
- **Implementation Quality**: Professional-grade code with comprehensive validation

### **Technical Specifications**
| Flux Scheme | Order | Shock Resolution | Contact Preservation | Computational Cost | Production Ready |
|-------------|-------|------------------|---------------------|-------------------|------------------|
| **Van Leer** | 1st | 4-5 cells | Exact | ~50 FLOPS/face | ✅ |
| **Roe (1st)** | 1st | 3-4 cells | Excellent | ~70 FLOPS/face | ✅ |
| **Roe (2nd)** | 2nd | 2-3 cells | Excellent | ~120 FLOPS/face | ✅ |
| **AUSM** | 1st | 3-4 cells | Good | ~60 FLOPS/face | ✅ |

### **Mathematical Framework**
All flux schemes implement rigorous mathematical foundations:
- **Hyperbolic Conservation Laws**: Proper treatment of Euler/Navier-Stokes equations
- **Riemann Problem Solutions**: Exact or approximate Riemann solvers with physical consistency
- **Entropy Conditions**: Thermodynamically admissible solutions with entropy fixes
- **TVD Properties**: Total Variation Diminishing schemes preventing spurious oscillations
- **Characteristic-Based Upwinding**: Proper wave decomposition and upwind bias

For comprehensive technical details, see the documentation files in the `docs/` directory.

## �📊 GPU Performance

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

#### **Core Solver Components**
- `Solver`: Main solver class with time integration and convergence control
- `Cell`: Computational cell representation with conservative/primitive variables
- `Face`: Face/interface handling with geometric properties and connectivity
- `Boundary_Conditions`: Comprehensive boundary condition implementations

#### **Flux Computation Framework** 🚀
- `Van_Leer()`: Van Leer flux vector splitting implementation in `src/Van_Leer.cpp`
- `ROE()`: Enhanced first-order Roe scheme with entropy fix in `src/Roe_Scheme.cpp`
- `ROE_2O()`: Second-order Roe scheme with TVD limiting in `src/Roe_Scheme.cpp`
- `AUSM()`: AUSM flux computation scheme in `src/Ausm_Flux.cpp`
- `Second_Order_Limiter()`: TVD slope limiters for high-resolution methods
- `Calculate_Primitive_Variables()`: Conservative to primitive variable conversion

#### **Advanced Features**
- **Entropy Fix**: Sonic point regularization preventing expansion shocks
- **Error Handling**: Comprehensive validation and graceful failure recovery
- **State Validation**: Physical consistency checks for density, pressure, temperature
- **Boundary Treatment**: Robust handling of wall, inlet, outlet, and far-field conditions

## 🔬 Development

### Code Structure
- **Object-Oriented Design**: Modular, maintainable code
- **CUDA Integration**: Seamless CPU-GPU data transfer
- **Memory Management**: Efficient memory usage patterns
- **Error Handling**: Comprehensive error checking

### Adding New Features

#### **Flux Scheme Development** 🔧
1. **New Flux Schemes**: 
   - Create dedicated source file in `src/` (e.g., `src/NewFlux_Scheme.cpp`)
   - Follow the established pattern with function signature: `void NEW_FLUX(int Cell_No, int N_Cell_No, int Face_No)`
   - Implement proper error checking and boundary condition handling
   - Add comprehensive documentation following existing templates

2. **Flux Scheme Integration**:
   - Update `Evaluate_Cell_Net_Flux_1O()` or `Evaluate_Cell_Net_Flux_2O()` functions
   - Add new dissipation type option in solver configuration
   - Include scheme in test case validation framework

#### **General Development Guidelines**
3. **New Numerical Schemes**: Add to `src/Numerical_Method.cpp` with proper mathematical documentation
4. **New Boundary Conditions**: Extend `src/Boundary_Conditions.cpp` with physical justification
5. **New CUDA Kernels**: Add to appropriate `CUDA_KERNELS/*.cu` file with performance optimization
6. **New Test Cases**: Create in `Test_Cases/` directory with analytical solutions for validation

#### **Documentation Standards** 📚
- **Mathematical Framework**: Include governing equations and derivations
- **Implementation Details**: Algorithm steps and computational complexity
- **Validation Results**: Comparison with analytical/reference solutions
- **Performance Analysis**: Computational cost and accuracy assessment
- **Usage Examples**: Clear examples with expected results

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
