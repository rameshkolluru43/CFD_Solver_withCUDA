# Testing and Validation Plan for New CUDA Kernels
**Date:** January 15, 2026  
**Status:** Ready for Testing on CUDA-enabled System

---

## System Requirements

### Required Hardware
- NVIDIA GPU with Compute Capability ≥ 6.0 (Pascal or newer)
- Recommended: RTX 2060 or better
- Minimum: GTX 1060 or better

### Required Software
- CUDA Toolkit ≥ 11.0
- CMake ≥ 3.16
- GCC/Clang with C++17 support
- VTK 9.x (for visualization)

### Verification
```bash
# Check CUDA availability
nvcc --version

# Check GPU
nvidia-smi

# Check CMake
cmake --version
```

---

## Compilation Steps

### Step 1: Clean Build
```bash
cd /Users/rameshkolluru/My_Research/CFD_Solver_withCUDA
rm -rf build/*
cd build
```

### Step 2: Configure with CMake
```bash
cmake .. -DCMAKE_BUILD_TYPE=Release
```

**Expected output:**
```
-- The CXX compiler identification is GNU/Clang
-- The CUDA compiler identification is NVIDIA
-- Detecting CUDA compile features - done
-- Configuring done
-- Generating done
```

**If errors occur:**
- Set CUDA path: `export CUDA_PATH=/usr/local/cuda`
- Set compiler: `cmake .. -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda`

### Step 3: Compile
```bash
make -j8 CFD_solver_gpu
```

**Expected output:**
```
[  5%] Building CUDA object CMakeFiles/...
[ 10%] Building CUDA object CMakeFiles/...
...
[100%] Linking CUDA executable CFD_solver_gpu
[100%] Built target CFD_solver_gpu
```

**Check for warnings:**
- Register usage warnings (acceptable if < 64 registers)
- Shared memory warnings (should not appear)
- Uncoalesced access warnings (investigate if present)

### Step 4: Verify Executable
```bash
ls -lh CFD_solver_gpu
./CFD_solver_gpu --help
```

---

## Unit Testing

### Test 1: Compilation Verification

**Purpose:** Ensure all kernels compile without errors

```bash
# Check if all object files were created
ls -lh CMakeFiles/CFD_solver_gpu.dir/CUDA_KERNELS/

# Should see:
# - Roe_Flux_Cuda_Kernels.cu.o
# - HLLC_LLF_Flux_Cuda_Kernels.cu.o
# - MUSCL_WENO_Reconstruction_Cuda_Kernels.cu.o
# - Advanced_Flux_Schemes_Cuda_Wrappers.cu.o
# - Reconstruction_Schemes_Cuda_Wrappers.cu.o
```

**Success criteria:** All .cu.o files present, no compilation errors

---

### Test 2: Sod Shock Tube (1D) - Roe Flux

**Purpose:** Validate Roe flux scheme against analytical solution

**Test code:** (Add to `Unit_Test_Codes/test_roe_sod.cpp`)

```cpp
#include "Advanced_Flux_Schemes_Cuda.h"
#include <cmath>
#include <iostream>

void test_roe_sod_shock() {
    // Parameters
    const int N = 1000;
    const double L = 1.0;
    const double dx = L / N;
    const double t_final = 0.2;
    const double gamma = 1.4;
    const double CFL = 0.5;
    
    // Initial conditions (Sod shock tube)
    double *U_cells = new double[N * 4];
    for (int i = 0; i < N; i++) {
        double x = i * dx;
        if (x < 0.5) {
            // Left state
            U_cells[i*4 + 0] = 1.0;        // rho
            U_cells[i*4 + 1] = 0.0;        // rho*u
            U_cells[i*4 + 2] = 0.0;        // rho*v
            U_cells[i*4 + 3] = 2.5;        // E
        } else {
            // Right state
            U_cells[i*4 + 0] = 0.125;      // rho
            U_cells[i*4 + 1] = 0.0;        // rho*u
            U_cells[i*4 + 2] = 0.0;        // rho*v
            U_cells[i*4 + 3] = 0.25;       // E
        }
    }
    
    // Allocate GPU memory
    double *d_U_cells, *d_P_cells, *d_flux;
    cudaMalloc(&d_U_cells, N * 4 * sizeof(double));
    cudaMalloc(&d_P_cells, N * 4 * sizeof(double));
    cudaMalloc(&d_flux, N * 4 * sizeof(double));
    
    // Copy to GPU
    cudaMemcpy(d_U_cells, U_cells, N*4*sizeof(double), cudaMemcpyHostToDevice);
    
    // Time integration loop
    double t = 0.0;
    int iter = 0;
    while (t < t_final) {
        // Compute primitives
        compute_primitives_kernel<<<(N+255)/256, 256>>>(d_U_cells, d_P_cells, N, gamma);
        
        // Compute Roe flux
        cudaError_t err = launch_roe_flux(d_U_cells, d_P_cells, 
                                          /*face data*/, d_flux, N, gamma);
        if (err != cudaSuccess) {
            std::cerr << "Roe flux failed: " << cudaGetErrorString(err) << std::endl;
            return;
        }
        
        // Time integration (RK4 or Euler)
        // ... update U_cells
        
        t += dt;
        iter++;
    }
    
    // Copy back to CPU
    cudaMemcpy(U_cells, d_U_cells, N*4*sizeof(double), cudaMemcpyDeviceToHost);
    
    // Validate against analytical solution
    // Expected at t=0.2:
    // - Shock at x ~ 0.85
    // - Contact at x ~ 0.69
    // - Rarefaction head at x ~ 0.26
    
    bool shock_ok = (U_cells[850*4] > 0.3);  // Shock density jump
    bool contact_ok = (U_cells[690*4] > 0.25 && U_cells[690*4] < 0.4);
    bool rarefaction_ok = (U_cells[260*4] < 0.9);
    
    std::cout << "Sod Shock Tube Test (Roe):" << std::endl;
    std::cout << "  Shock position: " << (shock_ok ? "PASS" : "FAIL") << std::endl;
    std::cout << "  Contact position: " << (contact_ok ? "PASS" : "FAIL") << std::endl;
    std::cout << "  Rarefaction: " << (rarefaction_ok ? "PASS" : "FAIL") << std::endl;
    
    // Cleanup
    delete[] U_cells;
    cudaFree(d_U_cells);
    cudaFree(d_P_cells);
    cudaFree(d_flux);
}
```

**Success criteria:**
- ✅ Shock at x ≈ 0.85 (± 0.05)
- ✅ Contact at x ≈ 0.69 (± 0.03)
- ✅ Rarefaction head at x ≈ 0.26 (± 0.03)
- ✅ No oscillations near discontinuities
- ✅ Pressure/density/velocity within 5% of analytical

---

### Test 3: Contact Discontinuity - HLLC Flux

**Purpose:** Verify HLLC resolves contacts exactly

```cpp
void test_hllc_contact() {
    // Initial: density jump only
    // Left:  rho=1.0, u=0.5, p=1.0
    // Right: rho=0.5, u=0.5, p=1.0
    
    // Run with HLLC
    launch_hllc_flux(...);
    
    // Validate: Contact should remain sharp
    // Check: No smearing of density profile
}
```

**Success criteria:**
- ✅ Contact discontinuity remains sharp (< 2 cells smearing)
- ✅ Velocity constant across contact
- ✅ Pressure constant across contact

---

### Test 4: Smooth Flow - WENO5 Order Verification

**Purpose:** Verify 5th-order convergence in smooth regions

```cpp
void test_weno5_convergence() {
    // Initial: u = sin(2*pi*x)
    // Exact solution: u(x,t) = sin(2*pi*(x-t))
    
    // Test with different resolutions: 50, 100, 200, 400
    double errors[4];
    int N[] = {50, 100, 200, 400};
    
    for (int i = 0; i < 4; i++) {
        // Run simulation
        launch_weno5_reconstruction(...);
        
        // Compute L2 error
        errors[i] = compute_L2_error(...);
    }
    
    // Compute convergence rate
    for (int i = 1; i < 4; i++) {
        double rate = log(errors[i-1]/errors[i]) / log(2.0);
        std::cout << "Convergence rate: " << rate << std::endl;
        // Should see rate ≈ 5.0 for WENO5
    }
}
```

**Success criteria:**
- ✅ Convergence rate ≈ 5.0 (± 0.5) for WENO5
- ✅ Convergence rate ≈ 2.0 (± 0.2) for MUSCL
- ✅ No oscillations in smooth regions

---

### Test 5: 2D Bow Shock - Robustness Test

**Purpose:** Test all schemes on complex 2D problem

```cpp
void test_bow_shock() {
    // Supersonic flow over cylinder
    // Mach = 2.0, gamma = 1.4
    
    // Test each scheme
    const char* schemes[] = {"Roe", "HLLC", "LLF"};
    
    for (const char* scheme : schemes) {
        std::cout << "Testing " << scheme << " scheme..." << std::endl;
        
        if (strcmp(scheme, "Roe") == 0)
            launch_roe_flux(...);
        else if (strcmp(scheme, "HLLC") == 0)
            launch_hllc_flux(...);
        else
            launch_llf_flux(...);
        
        // Check convergence
        bool converged = run_simulation();
        std::cout << scheme << ": " << (converged ? "PASS" : "FAIL") << std::endl;
    }
}
```

**Success criteria:**
- ✅ All schemes converge
- ✅ Shock standoff distance within 10% of expected
- ✅ No carbuncle instability
- ✅ Roe and HLLC produce similar results

---

### Test 6: GPU vs CPU Validation

**Purpose:** Ensure GPU results match CPU reference

```cpp
void test_gpu_vs_cpu() {
    // Run same problem on CPU and GPU
    
    // CPU: Use existing implementations
    compute_roe_flux_cpu(...);  // From src/Roe_Scheme.cpp
    
    // GPU: Use new kernels
    launch_roe_flux(...);
    
    // Compare results
    double max_diff = 0.0;
    for (int i = 0; i < N*4; i++) {
        double diff = fabs(cpu_result[i] - gpu_result[i]);
        max_diff = std::max(max_diff, diff);
    }
    
    std::cout << "Max difference (CPU vs GPU): " << max_diff << std::endl;
    // Should be < 1e-10 (machine precision)
}
```

**Success criteria:**
- ✅ Maximum difference < 1e-10
- ✅ RMS difference < 1e-12
- ✅ Identical results for all test cases

---

## Performance Benchmarking

### Benchmark 1: Flux Computation Speed

```cpp
void benchmark_flux_schemes() {
    const int N = 100000;  // Large grid
    const int num_trials = 100;
    
    // Benchmark Roe
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < num_trials; i++) {
        launch_roe_flux(...);
    }
    cudaDeviceSynchronize();
    auto end = std::chrono::high_resolution_clock::now();
    double time_roe = std::chrono::duration<double>(end - start).count() / num_trials;
    
    // Repeat for HLLC, LLF
    // ...
    
    // CPU reference
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < num_trials; i++) {
        compute_roe_flux_cpu(...);
    }
    end = std::chrono::high_resolution_clock::now();
    double time_cpu = std::chrono::duration<double>(end - start).count() / num_trials;
    
    // Print speedups
    std::cout << "Roe GPU time: " << time_roe*1000 << " ms" << std::endl;
    std::cout << "Roe CPU time: " << time_cpu*1000 << " ms" << std::endl;
    std::cout << "Speedup: " << time_cpu/time_roe << "x" << std::endl;
}
```

**Expected results:**
| Scheme | GPU Time (ms) | CPU Time (ms) | Speedup |
|--------|--------------|---------------|---------|
| Roe    | ~0.5         | ~20           | 40x     |
| HLLC   | ~0.3         | ~15           | 50x     |
| LLF    | ~0.2         | ~10           | 50x     |
| MUSCL  | ~0.4         | ~15           | 37x     |
| WENO5  | ~1.5         | ~40           | 27x     |

---

### Benchmark 2: Block Size Optimization

```cpp
void find_optimal_block_size() {
    int block_sizes[] = {32, 64, 128, 256, 512, 1024};
    
    for (int bs : block_sizes) {
        auto start = std::chrono::high_resolution_clock::now();
        launch_roe_flux(..., bs);
        cudaDeviceSynchronize();
        auto end = std::chrono::high_resolution_clock::now();
        
        double time = std::chrono::duration<double>(end - start).count();
        std::cout << "Block size " << bs << ": " << time*1000 << " ms" << std::endl;
    }
}
```

**Expected optimal block sizes:**
- Ampere (RTX 30xx): 256-512
- Turing (RTX 20xx): 256
- Pascal (GTX 10xx): 128-256

---

### Benchmark 3: Profiling with Nsight

```bash
# Profile Roe flux kernel
nsys profile --stats=true ./CFD_solver_gpu

# Output will show:
# - Kernel execution time
# - Memory transfer time
# - Occupancy
# - Register usage

# Check for:
# ✅ High occupancy (> 50%)
# ✅ Minimal __syncthreads() stalls
# ✅ Coalesced memory access
# ✅ No bank conflicts
```

---

## Integration Testing

### Full Solver Test

```bash
# Test with existing grid
./CFD_solver_gpu --input json_Files/cylinder_flow.json --flux-scheme Roe --reconstruction MUSCL

# Expected output:
# Iteration 100: Residual = 1.234e-03
# Iteration 200: Residual = 5.678e-05
# ...
# Converged in 500 iterations
```

**Success criteria:**
- ✅ Solver converges
- ✅ Results physically reasonable
- ✅ No NaN or Inf values
- ✅ Matches CPU results

---

## Debugging Guide

### Common Issues

**Issue 1: Illegal memory access**
```
Symptom: cuda error 77 (illegal memory access)
Cause: Invalid neighbor indices
Fix: Check face_neighbors and cell_neighbors arrays
Debug: Add bounds checking in kernel
```

**Issue 2: Solution diverges**
```
Symptom: NaN or Inf values
Cause: Numerical instability
Fix: 
  - Reduce CFL number
  - Use more robust scheme (LLF)
  - Enable entropy fix
```

**Issue 3: Wrong results**
```
Symptom: Results differ from CPU
Cause: 
  - Wrong data layout
  - Missing synchronization
  - Incorrect indexing
Fix: Add validation checks, print intermediate values
```

**Issue 4: Poor performance**
```
Symptom: Speedup < 10x
Cause:
  - Uncoalesced memory access
  - Too many registers
  - Low occupancy
Fix: Profile with nsys, optimize memory access
```

---

## Validation Checklist

Before marking implementation as production-ready:

### Correctness
- [ ] All unit tests pass
- [ ] GPU matches CPU results (< 1e-10 difference)
- [ ] Sod shock tube validation passed
- [ ] 2D test cases validated
- [ ] No NaN/Inf in any test

### Performance
- [ ] Speedup ≥ 20x for Roe flux
- [ ] Speedup ≥ 30x for HLLC flux
- [ ] Speedup ≥ 15x for WENO5
- [ ] Optimal block size determined
- [ ] Profiling completed

### Robustness
- [ ] No crashes on edge cases
- [ ] Handles boundary conditions correctly
- [ ] Works with different grid sizes
- [ ] Convergence for complex problems
- [ ] No memory leaks

### Documentation
- [ ] Integration guide followed
- [ ] All tests documented
- [ ] Performance metrics recorded
- [ ] Known limitations documented

---

## Final Steps

After all tests pass:

1. **Update Documentation**
   ```bash
   # Update the roadmap
   vi CUDA_IMPLEMENTATION_STATUS_AND_ROADMAP.md
   # Change status from 92% to 95%
   ```

2. **Commit Changes**
   ```bash
   git add CUDA_KERNELS/*.cu CUDA_KERNELS/*.h
   git add CMakeLists.txt
   git commit -m "Add Priority 1 & 2 CUDA kernels: Roe, HLLC, LLF, MUSCL, WENO5"
   git push
   ```

3. **Create Release Tag**
   ```bash
   git tag -a v2.0-cuda-complete -m "Complete Roe, HLLC, LLF flux schemes and MUSCL/WENO5 reconstruction"
   git push origin v2.0-cuda-complete
   ```

4. **Update README**
   - Add performance benchmarks
   - List supported schemes
   - Add usage examples

---

## Next: Priority 3 Implementation

Once validation is complete, proceed to Priority 3:
1. 3D extension of all schemes
2. Multi-GPU support
3. Advanced optimizations

---

**Status:** Ready for testing on CUDA-enabled system  
**Estimated Testing Time:** 2-4 hours  
**Estimated Validation Time:** 1-2 days
