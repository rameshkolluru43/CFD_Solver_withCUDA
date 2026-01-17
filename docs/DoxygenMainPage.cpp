/**
 * @mainpage CFD Solver Suite - Advanced Compressible/Incompressible Flow Solver with Turbulence Models
 *
 * @section intro_sec Introduction
 *
 * Welcome to the comprehensive technical documentation for the **Advanced CFD Solver Suite with CUDA Acceleration**.
 * This is a high-performance computational fluid dynamics solver featuring **production-ready implementations**
 * of advanced flux computation schemes, incompressible flow solver with SIMPLE algorithm, RANS turbulence models,
 * and GPU acceleration capabilities for both compressible and incompressible flows.
 *
 * @section version_info Version Information
 *
 * - **Version**: v3.0 - Incompressible Solver & Turbulence Models
 * - **Release Date**: September 2025
 * - **Language**: C++ with CUDA support
 * - **License**: GNU General Public License v3.0
 *
 * @section key_features Key Features
 *
 * ### 🚀 Compressible Flow: Advanced Flux Computation Schemes
 *
 * #### **Van Leer Flux Vector Splitting**
 * - Complete implementation with Mach number-based splitting
 * - Exact contact discontinuity preservation
 * - Excellent performance across all flow regimes
 * - ~50 floating point operations per face
 * - **Implementation**: Van_Leer() function in src/Van_Leer.cpp
 *
 * #### **Enhanced Roe Approximate Riemann Solver**
 * - **First-Order**: Comprehensive error checking, entropy fix, boundary handling
 * - **Second-Order**: TVD slope limiting with high-resolution shock capturing
 * - Mathematical rigor with complete eigenvalue-eigenvector decomposition
 * - Production-ready reliability with industrial-grade error handling
 * - **Implementation**: ROE() and ROE_2O() functions in src/Roe_Scheme.cpp
 *
 * #### **AUSM Flux Scheme**
 * - All-speed capability from incompressible to hypersonic regimes
 * - Advanced upwind splitting for momentum and energy equations
 * - Robust pressure and acoustic wave treatment
 * - **Implementation**: AUSM() function in src/Ausm_Flux.cpp
 *
 * ### 🌊 Incompressible Flow: SIMPLE Algorithm (NEW v3.0)
 *
 * #### **Cell-Centered Staggered Grid Finite Volume Method**
 * - **SIMPLE Algorithm**: Semi-Implicit Method for Pressure-Linked Equations
 * - **Pressure-Velocity Coupling**: Rhie-Chow interpolation for stability
 * - **Iterative Solvers**: BiCGSTAB and GMRES with preconditioning
 * - **Boundary Conditions**: Inlet velocity, outlet pressure, wall, symmetry
 * - **Convergence Control**: Adaptive relaxation and residual monitoring
 * - **Implementation**: IncompressibleSolver class in Incompressible_Solver/
 *
 * #### **Key Incompressible Features**
 * - Momentum equation discretization with central differencing
 * - Pressure correction equation with proper boundary conditions
 * - Under-relaxation for pressure and velocity fields
 * - Comprehensive boundary condition implementation
 * - Standalone and integrated operation modes
 *
 * ### 🌪️ RANS Turbulence Models (NEW v3.0)
 *
 * #### **K-epsilon Turbulence Model**
 * - **Standard Model**: Launder-Spalding formulation
 * - **Turbulent Viscosity**: Proper calculation with realizability constraints
 * - **Wall Functions**: Enhanced wall treatment for coarse grids
 * - **Source Terms**: Production and dissipation with proper limiting
 * - **Implementation**: K_Epsilon_Model class in src/K_Epsilon_Model.cpp
 *
 * #### **K-omega Turbulence Models**
 * - **Wilcox K-omega**: Original formulation with low-Re corrections
 * - **SST K-omega**: Menter's Shear Stress Transport model
 * - **Automatic Wall Distance**: Efficient wall distance calculation
 * - **Blending Functions**: Smooth transition between inner/outer regions
 * - **Cross-diffusion Terms**: Proper SST model implementation
 * - **Implementation**: K_Omega_Model class in src/K_Omega_Model.cpp
 *
 * #### **Turbulence Integration Features**
 * - Seamless integration with both compressible and incompressible solvers
 * - RK4 time integration for turbulence transport equations
 * - Proper source term linearization for stability
 * - Comprehensive boundary condition treatment
 * - Production limiting for numerical stability
 *
 * ### 🔧 Computational Framework
 *
 * - **Flow Physics**: Full Euler/Navier-Stokes (compressible) and SIMPLE (incompressible)
 * - **Time Integration**: Explicit Runge-Kutta (RK4), TVD-RK3, implicit schemes
 * - **Spatial Discretization**: High-order MUSCL reconstruction, central differencing
 * - **Slope Limiters**: Van Leer, Minmod, Superbee with TVD properties
 * - **Matrix Solvers**: BiCGSTAB, GMRES, direct solvers for incompressible systems
 * - **GPU Acceleration**: CUDA kernels for all major computations
 * - **Error Handling**: Comprehensive validation and graceful recovery
 *
 * @section architecture_overview Architecture Overview
 *
 * ### Core Components
 *
 * | Component | Description | Key Files |
 * |-----------|-------------|-----------|
 * | **Compressible Flux Computation** | Advanced Riemann solvers | src/Van_Leer.cpp, src/Roe_Scheme.cpp, src/Ausm_Flux.cpp |
 * | **Incompressible Solver** | SIMPLE algorithm implementation | Incompressible_Solver/Incompressible_Solver.cpp |
 * | **Turbulence Models** | RANS turbulence modeling | src/K_Epsilon_Model.cpp, src/K_Omega_Model.cpp |
 * | **Time Integration** | RK4, TVD-RK3, implicit schemes | src/Time_Integration.cpp, src/Turbulence_Integration.cpp |
 * | **Boundary Conditions** | Wall, inlet, outlet, far-field | src/Boundary_Conditions.cpp |
 * | **CUDA Kernels** | GPU-accelerated computations | CUDA_KERNELS/*.cu |
 * | **Grid Processing** | Structured/unstructured mesh handling | src/Grid_Functions.cpp |
 * | **Matrix Solvers** | Iterative linear system solvers | src/Matrix_Solvers.cpp |
 * | **I/O System** | VTK output and JSON configuration | src/IO_Functions.cpp |
 *
 * @section mathematical_framework Mathematical Framework
 *
 * ### Governing Equations
 *
 * The solver implements the compressible Navier-Stokes equations in conservative form:
 *
 * ∂U/∂t + ∇·F(U) = ∇·Fᵥ(U,∇U)
 *
 * Where:
 * - **U**: Conservative variable vector [ρ, ρu, ρv, ρE]ᵀ
 * - **F(U)**: Inviscid flux tensor (Euler fluxes)
 * - **Fᵥ(U,∇U)**: Viscous flux tensor (Navier-Stokes terms)
 *
 * ### Flux Computation Philosophy
 *
 * All flux schemes implement:
 * - **Hyperbolic Conservation Laws**: Proper treatment of wave propagation
 * - **Riemann Problem Solutions**: Exact or approximate with physical consistency
 * - **Entropy Conditions**: Thermodynamically admissible solutions
 * - **TVD Properties**: Total Variation Diminishing for shock capturing
 * - **Characteristic-Based Upwinding**: Proper wave decomposition
 *
 * @section performance_benchmarks Performance Benchmarks
 *
 * ### Flux Scheme Comparison
 *
 * | Scheme | Order | Shock Resolution | Contact Preservation | FLOPS/Face | Production Ready |
 * |--------|-------|------------------|---------------------|------------|------------------|
 * | Van Leer | 1st | 4-5 cells | Exact | ~50 | ✅ |
 * | Roe (1st) | 1st | 3-4 cells | Excellent | ~70 | ✅ |
 * | Roe (2nd) | 2nd | 2-3 cells | Excellent | ~120 | ✅ |
 * | AUSM | 1st | 3-4 cells | Good | ~60 | ✅ |
 *
 * ### GPU Performance
 * - **Speedup**: 10-50x over CPU implementation depending on problem size
 * - **Memory Bandwidth**: Optimized coalesced memory access patterns
 * - **Architecture Support**: CUDA compute capability 6.0-9.0 (Pascal to Hopper)
 *
 * @section validation_framework Validation Framework
 *
 * ### Standard Test Cases
 *
 * #### Compressible Flow Validation
 * - **Sod Shock Tube**: 1D Riemann problem validation
 * - **Lax Problem**: Strong shock and rarefaction interaction
 * - **Woodward-Colella**: Complex shock-shock interactions
 * - **Flow Over Cylinder**: 2D viscous boundary layer flows
 * - **Double Mach Reflection**: Shock-boundary interactions
 *
 * #### Incompressible Flow Validation
 * - **Lid-Driven Cavity**: Classical benchmark at various Reynolds numbers
 * - **Channel Flow**: Fully developed turbulent channel flow
 * - **Backward-Facing Step**: Separated flow with reattachment
 * - **Flow Over Cylinder**: Low-speed viscous flow validation
 *
 * #### Turbulence Model Validation
 * - **Flat Plate Boundary Layer**: Zero pressure gradient validation
 * - **Channel Flow**: Turbulent channel with wall functions
 * - **Backward-Facing Step**: Separated flow with turbulence
 * - **Mixing Layer**: Free shear layer turbulence
 *
 * ### Validation Metrics
 * - **L₁, L₂, L∞ Error Norms**: Quantitative accuracy assessment
 * - **Convergence Rates**: Spatial and temporal convergence analysis
 * - **Conservation Properties**: Mass, momentum, energy conservation
 * - **Entropy Production**: Physical admissibility verification
 *
 * @section usage_examples Usage Examples
 *
 * ### Compressible Flow Setup
 * @code{.cpp}
 * // Initialize solver with Van Leer flux scheme
 * Dissipation_Type = 1;  // Van Leer
 *
 * // Configure for second-order accuracy
 * Is_Second_Order = true;
 *
 * // Enable turbulence modeling
 * Enable_Turbulence = true;
 * Turbulence_Model = 1;  // K-epsilon model
 *
 * // Run simulation
 * for (int iter = 0; iter < Max_Iterations; iter++) {
 *     Evaluate_Cell_Net_Flux_2O();  // High-resolution flux evaluation
 *     UpdateTurbulenceVariables_RK4(dt);  // Turbulence time integration
 *     UpdateConservativeVariables_RK4(dt);  // Flow field time integration
 *
 *     if (Check_Convergence()) break;
 * }
 * @endcode
 *
 * ### Incompressible Flow Setup
 * @code{.cpp}
 * // Initialize SIMPLE solver
 * IncompressibleSolver solver;
 * solver.LoadConfiguration("incompressible_config.json");
 *
 * // Set boundary conditions
 * solver.SetInletVelocity(inlet_faces, velocity_profile);
 * solver.SetOutletPressure(outlet_faces, reference_pressure);
 * solver.SetWallBoundary(wall_faces);
 *
 * // Run SIMPLE iterations
 * for (int iter = 0; iter < max_iterations; iter++) {
 *     solver.SolveMomentumEquations();
 *     solver.SolvePressureCorrection();
 *     solver.UpdateVelocityField();
 *     solver.UpdatePressureField();
 *
 *     if (solver.CheckConvergence()) break;
 * }
 * @endcode
 *
 * ### Flux Scheme Selection
 * @code{.cpp}
 * // Available compressible flux schemes
 * Dissipation_Type = 1;  // Van Leer flux vector splitting
 * Dissipation_Type = 2;  // LLF (Local Lax-Friedrichs)
 * Dissipation_Type = 3;  // Enhanced Roe scheme (1st/2nd order)
 * Dissipation_Type = 4;  // AUSM flux scheme
 * @endcode
 *
 * ### Turbulence Model Selection
 * @code{.cpp}
 * // Available turbulence models
 * Turbulence_Model = 0;  // No turbulence (laminar)
 * Turbulence_Model = 1;  // K-epsilon model
 * Turbulence_Model = 2;  // K-omega (Wilcox)
 * Turbulence_Model = 3;  // SST K-omega
 * @endcode
 *
 * @section documentation_structure Documentation Structure
 *
 * ### API Reference
 * - **Classes**: Core solver components and data structures
 * - **Functions**: Detailed function documentation with parameters
 * - **Files**: Source code organization and module descriptions
 * - **Namespaces**: Logical code organization
 *
 * ### Technical Guides
 * - **docs/Van_Leer_Flux_Implementation.md**: Van Leer scheme mathematical framework
 * - **docs/ROE_2O_Implementation.md**: Second-order Roe implementation details
 * - **docs/Enhanced_ROE_First_Order.md**: Enhanced first-order Roe documentation
 * - **docs/AUSM_Flux_Implementation.md**: AUSM flux scheme technical guide
 *
 * ### Completion Summaries
 * - Implementation achievement summaries with before/after comparisons
 * - Performance analysis and validation results
 * - Usage examples and integration guidelines
 *
 * @section getting_started Getting Started
 *
 * ### Prerequisites
 * - **CMake** (≥3.16): Build system
 * - **CUDA Toolkit** (≥11.0): GPU computation
 * - **C++ Compiler**: C++17 compatible (GCC/Clang)
 * - **JsonCpp**: Configuration file parsing
 * - **VTK** (≥9.4): Visualization output
 *
 * ### Quick Build
 * @code{.bash}
 * mkdir build && cd build
 * cmake ..
 * make -j$(nproc)
 * @endcode
 *
 * ### Running Your First Simulation
 * @code{.bash}
 * # CPU version
 * ./CFD_solver ../json_Files/Solver_Config.json
 *
 * # GPU version (if CUDA available)
 * ./CFD_solver_gpu ../json_Files/Solver_Config.json
 * @endcode
 *
 * @section contributing Contributing
 *
 * ### Development Standards
 * - **Mathematical Documentation**: Include governing equations and derivations
 * - **Code Quality**: Comprehensive error checking and validation
 * - **Performance**: Optimize for both accuracy and computational efficiency
 * - **Testing**: Validate against analytical/reference solutions
 *
 * ### Flux Scheme Development
 * 1. Create dedicated source file following naming convention
 * 2. Implement standard function signature with proper error handling
 * 3. Add comprehensive mathematical documentation
 * 4. Include validation test cases and performance analysis
 * 5. Update configuration system and user documentation
 *
 * @section references References
 *
 * ### Key Literature
 * - **Toro, E.F.**: "Riemann Solvers and Numerical Methods for Fluid Dynamics"
 * - **LeVeque, R.J.**: "Finite Volume Methods for Hyperbolic Problems"
 * - **Hirsch, C.**: "Numerical Computation of Internal and External Flows"
 * - **Blazek, J.**: "Computational Fluid Dynamics: Principles and Applications"
 *
 * ### Flux Scheme References
 * - **Van Leer, B.**: "Flux Vector Splitting for the Euler Equations" (1982)
 * - **Roe, P.L.**: "Approximate Riemann Solvers" (1981)
 * - **Liou, M.S.**: "A Sequel to AUSM: AUSM+" (1996)
 * - **Harten, A.**: "High Resolution Schemes for Hyperbolic Conservation Laws" (1983)
 *
 * @section contact_info Contact Information
 *
 * **Developer**: Ramesh Kolluru
 * **Email**: rameshkolluru43@gmail.com
 * **Project**: https://github.com/rameshkolluru43/CFD_Solver_withCUDA
 * **License**: GPL v3.0
 *
 * ---
 *
 * *This solver represents state-of-the-art computational fluid dynamics with production-ready
 * flux computation schemes, comprehensive error handling, and GPU acceleration for
 * high-performance scientific computing applications.*
 */