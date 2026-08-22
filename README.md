# CFD Solver with CUDA Kernels

A high-performance Computational Fluid Dynamics (CFD) solver featuring GPU acceleration through CUDA kernels. The solver supports both Euler and Navier-Stokes equations with various numerical schemes and boundary conditions for compressible flow simulations.

**Start here:** [LAYOUT.md](LAYOUT.md) · [docs/overview.md](docs/overview.md) · [docs/README.md](docs/README.md)

### Recent

- **Folder layout:** code in `solvers/`, JSON/meshes in `input/`, VTK/plots in `solutions/` (see [LAYOUT.md](LAYOUT.md)).
- **3D path (`solvers/src3d/`):** Cartesian Euler/NS, WENO-Z, WALE LES, SWBLI plate, Mach-6 open cavity (90° block mesh). Binaries `CFD_solver_3d` / `CFD_solver_3d_gpu`.
- **Half-cylinder viscous corner fix (July 2026):** LBRT face normals. Plots: `solutions/plots/half_cylinder_inviscid_viscous/`. Details: [docs/overview.md](docs/overview.md), [docs/RELEASE_NOTES.md](docs/RELEASE_NOTES.md).

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
- **CUDA Main-Loop Flux Path**: `CFD_solver_gpu` includes a CUDA 1st-order inviscid net-flux dispatcher for **LLF**, **MOVERS**, **RICCA**, and **MOVERS_NWSC** (`Dissipation_Type` 1, 2, 4, 5)
- **CUDA Kernels**: Additional kernels are available for Roe/HLLC/LLF flux experiments, reconstruction, gradients, viscous flux prototypes, matrix assembly, and time-integration utilities
- **Memory Management**: Efficient host-device memory transfers
- **Parallel Execution**: CUDA architecture is configured in `CMakeLists.txt` (`CUDA_ARCHITECTURES`; currently tuned for compute 7.5 and adjustable for your GPU)
- **Iterative Solvers**: GPU-accelerated linear algebra operations
- **Metal (experimental, macOS / Apple Silicon)**: Optional `solvers/Metal_Kernels/` bridge and `.metal` shaders — standalone `Makefile` in that directory

Current CUDA main-loop limits: 2nd-order MUSCL flux, WENO flux, viscous flux, Runge-Kutta state updates, implicit solving, and AMR split/merge still execute through CPU-side loop logic.

### CPU Parallelism
- **OpenMP**: Multi-threaded host loops (`-fopenmp` / `libgomp` on Linux; Homebrew `libomp` on macOS)
- **MPI**: Optional `CFD_solver_mpi` — replicated mesh, rank-owned cell ranges, allgather of conservatives (`mpirun --oversubscribe -np 2` on small machines)

### Solver Capabilities

#### **Flow Physics** 🌊
- **Compressible Flow**: Full Euler and Navier-Stokes equations with thermodynamic consistency
- **Flow Regimes**: Subsonic, transonic, supersonic, and hypersonic flow support
- **Shock Capturing**: Advanced shock-capturing with entropy-satisfying schemes
- **Contact Preservation**: Exact contact discontinuity preservation with specialized flux methods

#### **Computational Framework** 🔧
- **Boundary Conditions**: Far-field, wall, inlet, outlet, symmetry with robust treatment
- **Grid Support**: Unstructured grids via **VTK/GMSH**; **mixed triangle and quadrilateral** 2D meshes with per-face connectivity
- **Dynamic Meshing**: **Gradient-based adaptive mesh refinement (AMR)** — tag cells by density/pressure gradients; configurable threshold and period
- **Test Cases**: Comprehensive validation suite including shock tubes, cylinder flows, and complex geometries
- **Output Formats**: VTK for ParaView visualization with field variable export
- **Error Handling**: Production-ready error checking and graceful failure recovery

## 📁 Project Structure

See **[LAYOUT.md](LAYOUT.md)** for binaries, cwd rules, and example commands.

```
solvers/src, src3d, include, CUDA_KERNELS   # 2D + 3D + CUDA
input/json_Files, Grid_Files                # run JSON and 2D meshes
wrappers/scripts, python                    # helpers
tests/Test_Cases                            # 2D case drivers
solutions/2D_Euler_Solutions                # VTK + case plots
docs/                                       # guides
build/                                      # CMake binaries
```

## 🛠️ Dependencies

### 2D CPU (`CFD_solver`)
- **CMake** ≥ 3.16, **C++17**
- **Boost** (regex), **JsonCpp**
- **OpenMP** (typical with GCC)
- **VTK** optional (library VTK I/O). Structured **TXT** meshes run without VTK.

### 2D MPI (`CFD_solver_mpi`)
- OpenMPI runtime (`mpirun`, `libmpi`) plus `mpi.h` (system `libopenmpi-dev` or a user prefix such as `$HOME/deps`)

### GPU (`CFD_solver_gpu`, `CFD_solver_3d_gpu`)
- **CUDA Toolkit** with `nvcc` on `PATH`. Architecture: `CFD_CUDA_ARCH` / `CFD_CUDA_ARCHITECTURES` (default **75**)

### 3D (`CFD_solver_3d`)
- Same C++/OpenMP stack; no VTK library required (writes legacy VTK text)

### Optional
- **Doxygen**, **GMSH**, **ParaView**

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
git clone https://github.com/rameshkolluru43/CFD_Solver_withCUDA.git
cd CFD_Solver_withCUDA

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

Typical binaries in `build/`: `CFD_solver`, `CFD_solver_mpi`, `CFD_solver_gpu`, `CFD_solver_3d`, `CFD_solver_3d_gpu` (each only if its dependencies were found).

## 🚀 Usage

### Running Simulations

3D configs use `VTK_Out` under `solutions/2D_Euler_Solutions/`. Run **3D from the repo root**. Run **2D from `build/`** so `../input/Grid_Files` resolves. Full table: [LAYOUT.md](LAYOUT.md).

#### 3D
```bash
./build/CFD_solver_3d input/json_Files/Run_3D_Sod_LLF_smoke.json
./build/CFD_solver_3d_gpu input/json_Files/Run_3D_Cavity_M6_LES_WALE.json
```

#### 2D CPU / OpenMP / MPI / GPU
```bash
cd build
./CFD_solver ../input/json_Files/Solver_Config.json
OMP_NUM_THREADS=4 ./CFD_solver ../input/json_Files/MPI_Smoke_Test.json
mpirun --oversubscribe -np 2 ./CFD_solver_mpi ../input/json_Files/MPI_Smoke_Test.json
./CFD_solver_gpu ../input/json_Files/Run_HalfCylinder_M6_WENO_wallfix_smoke50.json
```

`CFD_solver_gpu` follows the same solver configuration as the CPU executable. For inviscid first-order runs, `Dissipation_Type` selects the CUDA flux kernel when available:

| `Dissipation_Type` | Scheme | CUDA main-loop support |
|--------------------|--------|------------------------|
| 1 | LLF | Yes |
| 2 | MOVERS | Yes |
| 3 | Roe | CPU fallback in current main-loop dispatcher |
| 4 | RICCA | Yes |
| 5 | MOVERS_NWSC | Yes |

See [docs/BUILD_AND_RUN.md](docs/BUILD_AND_RUN.md) for detailed build and run instructions.

### Configuration

The solver uses JSON configuration files to define simulation parameters:

#### Main Configuration (`Solver_Config.json`)
```json
{
    "TestCase": {
        "Test_Case": 1,
        "Test_Case_Name": "Half_Cylinder",
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
6. **SWBLI**: 2D and 3D shock-wave / boundary-layer interaction
7. **3D ramp / cavity**: 15° ramp and Mach-6 open cavity (`solvers/src3d/`)
8. **Scramjet Inlet**: Hypersonic inlet flow

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

For detailed testing information, see [docs/TEST_SUMMARY.md](docs/TEST_SUMMARY.md).

## 📊 Flux Scheme Enhancements

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

## 📊 GPU Performance

The CUDA implementation currently accelerates selected solver paths and provides kernels for future full-GPU integration:
- **Main-loop inviscid flux**: CUDA 1st-order net-flux path for LLF, MOVERS, RICCA, and MOVERS_NWSC.
- **Host/device bridge**: Current bridge flattens the existing CPU-side `Cell`/state vectors, launches CUDA, and copies back `Cells_Net_Flux` and `del_t`.
- **Kernel library**: Roe, HLLC/LLF, MUSCL/WENO reconstruction, gradient, viscous, matrix assembly, and time integration kernels are available for integration work.
- **Remaining CPU-side paths**: 2nd-order MUSCL dispatch, WENO flux in the active loop, RK updates, viscous fluxes, AMR topology changes, and implicit solving.

### Supported CUDA Architectures
- Pascal (6.0, 6.1)
- Volta (7.0)
- Turing (7.5)
- Ampere (8.0, 8.6)
- Ada Lovelace (8.9)
- Hopper (9.0)

## 📖 Documentation

| Document | Description |
|----------|-------------|
| [docs/overview.md](docs/overview.md) | Project overview, CUDA/host status, half-cylinder M=6 |
| [LAYOUT.md](LAYOUT.md) | Folder map, binaries, cwd, example commands |
| [docs/README.md](docs/README.md) | Full documentation index |
| [docs/HALF_CYLINDER_VALIDATION.md](docs/HALF_CYLINDER_VALIDATION.md) | Canonical M=6 configs, expected Pmax/Mmax, plots |
| [docs/MESH_AND_GRID.md](docs/MESH_AND_GRID.md) | Mesh formats, **LBRT** structured TXT, face normals, mixed tri/quad |
| [docs/ADAPTIVE_MESH_REFINEMENT.md](docs/ADAPTIVE_MESH_REFINEMENT.md) | Gradient-based AMR: indicator, tagging, configuration |
| [docs/CONFIGURATION.md](docs/CONFIGURATION.md) | JSON configuration reference (solver, viscous, AMR) |
| [docs/BUILD_AND_RUN.md](docs/BUILD_AND_RUN.md) | Build requirements, CMake options, running CPU/GPU |
| [docs/RELEASE_NOTES.md](docs/RELEASE_NOTES.md) | Changelog and recent feature summary |
| [docs/MAIN_LOOP_STATUS.md](docs/MAIN_LOOP_STATUS.md) | Main loop: host vs CUDA, WENO, viscous, AMR |

### Generate Documentation
```bash
# From the repository root
./wrappers/scripts/update_docs.sh
# or:
doxygen docs/Doxyfile_Cleaned
```

HTML output: **`docs/doxygen/html/index.html`** (`OUTPUT_DIRECTORY` in `docs/Doxyfile_Cleaned`).

### Key Classes and Functions

#### **Core Solver Components**
- `Solver`: Main solver class with time integration and convergence control
- `Cell`: Computational cell representation with conservative/primitive variables
- `Face`: Face/interface handling with geometric properties and connectivity
- `Boundary_Conditions`: Comprehensive boundary condition implementations

#### **Flux Computation Framework** 🚀
- `Van_Leer()`: `solvers/src/Van_Leer.cpp`
- `ROE()` / `ROE_2O()`: `solvers/src/Roe_Scheme.cpp`
- `AUSM()`: `solvers/src/Ausm_Flux.cpp`
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
   - Create dedicated source file in `solvers/src/` (e.g. `solvers/src/NewFlux_Scheme.cpp`)
   - Follow the established pattern with function signature: `void NEW_FLUX(int Cell_No, int N_Cell_No, int Face_No)`
   - Implement proper error checking and boundary condition handling
   - Add comprehensive documentation following existing templates

2. **Flux Scheme Integration**:
   - Update `Evaluate_Cell_Net_Flux_1O()` or `Evaluate_Cell_Net_Flux_2O()` functions
   - Add new dissipation type option in solver configuration
   - Include scheme in test case validation framework

#### **General Development Guidelines**
3. **New Numerical Schemes**: Add to `solvers/src/Numerical_Method.cpp` with proper mathematical documentation
4. **New Boundary Conditions**: Extend `solvers/src/Boundary_Conditions.cpp` with physical justification
5. **New CUDA Kernels**: Add to `solvers/CUDA_KERNELS/*.cu`
6. **Metal (macOS)**: `solvers/Metal_Kernels/` — `scripts/build_metallib.sh`
7. **New Test Cases**: `tests/Test_Cases/`

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
