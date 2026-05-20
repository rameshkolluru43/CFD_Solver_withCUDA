# Integration Guide for New CUDA Flux and Reconstruction Schemes
**Date:** January 15, 2026  
**Status:** Implementation Complete - Ready for Integration

## Overview
This guide describes how to integrate the newly implemented Priority 1 and Priority 2 CUDA kernels into the main CFD solver.

---

## New Files Created

### Kernel Implementation Files
1. **Roe_Flux_Cuda_Kernels.cu** (~600 lines)
   - `roe_flux_kernel` - 1st-order Roe flux
   - `roe_flux_2nd_order_kernel` - 2nd-order Roe with MUSCL

2. **HLLC_LLF_Flux_Cuda_Kernels.cu** (~400 lines)
   - `hllc_flux_kernel` - HLLC flux scheme
   - `llf_flux_kernel` - Local Lax-Friedrichs flux

3. **MUSCL_WENO_Reconstruction_Cuda_Kernels.cu** (~700 lines)
   - `muscl_reconstruction_kernel` - 2nd-order MUSCL
   - `weno5_reconstruction_kernel` - 5th-order WENO (component-wise)
   - `weno5_characteristic_reconstruction_kernel` - 5th-order WENO (characteristic)

### Header Files
4. **Advanced_Flux_Schemes_Cuda.h**
   - Declarations for Roe, HLLC, LLF kernels and wrappers

5. **Reconstruction_Schemes_Cuda.h**
   - Declarations for MUSCL and WENO kernels and wrappers

### Wrapper Files
6. **Advanced_Flux_Schemes_Cuda_Wrappers.cu**
   - Host wrapper functions for flux schemes

7. **Reconstruction_Schemes_Cuda_Wrappers.cu**
   - Host wrapper functions for reconstruction schemes

---

## Compilation Instructions

### 1. Update CMakeLists.txt

Add the new CUDA source files to your CMakeLists.txt:

```cmake
# Add new CUDA kernel files
set(CUDA_SOURCES
    # ... existing sources ...
    CUDA_KERNELS/Roe_Flux_Cuda_Kernels.cu
    CUDA_KERNELS/HLLC_LLF_Flux_Cuda_Kernels.cu
    CUDA_KERNELS/MUSCL_WENO_Reconstruction_Cuda_Kernels.cu
    CUDA_KERNELS/Advanced_Flux_Schemes_Cuda_Wrappers.cu
    CUDA_KERNELS/Reconstruction_Schemes_Cuda_Wrappers.cu
)
```

### 2. Compilation Test

```bash
cd /Users/rameshkolluru/My_Research/CFD_Solver_withCUDA/build
cmake ..
make -j8
```

---

## Integration Steps

### Step 1: Include Headers in Main Solver

In your main solver file (e.g., `Main/solver.cpp` or `CUDA_INTEGRATED_NUMERICAL_METHOD.cpp`):

```cpp
#include "Advanced_Flux_Schemes_Cuda.h"
#include "Reconstruction_Schemes_Cuda.h"
```

### Step 2: Add Runtime Selection for Flux Schemes

Add an enum for flux scheme selection:

```cpp
enum FluxScheme {
    AUSM_PLUS,      // Existing
    VAN_LEER,       // Existing
    CENTRAL,        // Existing
    ROE,            // NEW
    HLLC,           // NEW
    LLF             // NEW
};

FluxScheme active_flux_scheme = ROE;  // Default to Roe
```

### Step 3: Modify Flux Computation Loop

Replace or augment your existing flux computation with scheme selection:

```cpp
// In your main solver loop
switch(active_flux_scheme) {
    case ROE:
        if (use_second_order) {
            cudaError = launch_roe_flux_2nd_order(
                d_U_cells, d_P_cells, d_face_neighbors,
                d_face_normals, d_face_areas,
                d_limited_gradients, d_cell_centers, d_face_centers,
                d_flux_output, num_cells, gamma
            );
        } else {
            cudaError = launch_roe_flux(
                d_U_cells, d_P_cells, d_face_neighbors,
                d_face_normals, d_face_areas,
                d_flux_output, num_cells, gamma
            );
        }
        break;
        
    case HLLC:
        cudaError = launch_hllc_flux(
            d_U_cells, d_P_cells, d_face_neighbors,
            d_face_normals, d_face_areas,
            d_flux_output, num_cells, gamma
        );
        break;
        
    case LLF:
        cudaError = launch_llf_flux(
            d_U_cells, d_P_cells, d_face_neighbors,
            d_face_normals, d_face_areas,
            d_flux_output, num_cells, gamma
        );
        break;
        
    // ... other schemes ...
}

// Check for errors
if (cudaError != cudaSuccess) {
    fprintf(stderr, "Flux computation failed: %s\n", 
            cudaGetErrorString(cudaError));
    exit(1);
}
```

### Step 4: Add Runtime Selection for Reconstruction

```cpp
enum ReconstructionScheme {
    FIRST_ORDER,    // No reconstruction
    MUSCL,          // NEW
    WENO5,          // NEW
    WENO5_CHAR      // NEW (characteristic-based)
};

ReconstructionScheme active_reconstruction = MUSCL;
```

### Step 5: Integrate Reconstruction in Solver

Before flux computation (for 2nd-order or higher schemes):

```cpp
if (active_reconstruction == MUSCL) {
    // First compute gradients and limiters (existing code)
    launch_compute_gradients(...);
    launch_apply_limiters(...);
    
    // Then apply MUSCL reconstruction
    cudaError = launch_muscl_reconstruction(
        d_U_cells, d_cell_neighbors, d_cell_distances,
        d_limited_gradients, d_Q_left, d_Q_right,
        num_cells, kappa  // kappa = 0.5 for QUICK scheme
    );
}
else if (active_reconstruction == WENO5) {
    cudaError = launch_weno5_reconstruction(
        d_U_cells, d_cell_neighbors,
        d_Q_left, d_Q_right, num_cells
    );
}
else if (active_reconstruction == WENO5_CHAR) {
    // Requires pre-computed transformation matrices
    cudaError = launch_weno5_characteristic_reconstruction(
        d_U_cells, d_cell_neighbors,
        d_L_matrices, d_R_matrices,
        d_Q_left, d_Q_right, num_cells
    );
}
```

### Step 6: Configuration File Support

Add to your JSON configuration file:

```json
{
  "solver": {
    "flux_scheme": "Roe",
    "reconstruction": "MUSCL",
    "order": 2,
    "muscl_kappa": 0.5,
    "weno_epsilon": 1e-6
  }
}
```

---

## Testing and Validation

### Test Case 1: Sod Shock Tube (1D verification)

```cpp
void test_roe_sod_shock() {
    // Initial conditions:
    // Left:  rho=1.0, u=0.0, p=1.0
    // Right: rho=0.125, u=0.0, p=0.1
    // Diaphragm at x=0.5
    
    // Expected results at t=0.2:
    // - Shock at x ~ 0.85
    // - Contact at x ~ 0.69
    // - Rarefaction head at x ~ 0.26
    
    // Compare with analytical solution
}
```

### Test Case 2: Contact Discontinuity

```cpp
void test_hllc_contact() {
    // HLLC should resolve contacts exactly
    // Initial: jump in density only
    // Check that contact remains sharp
}
```

### Test Case 3: Smooth Flow (Order Verification)

```cpp
void test_weno5_smooth() {
    // Use smooth initial condition: sin(2πx)
    // Measure convergence rate
    // Should see ~5th order for WENO5
}
```

### Test Case 4: Bow Shock (2D test)

```cpp
void test_roe_bow_shock() {
    // Supersonic flow over cylinder
    // Verify shock capturing
    // Compare Roe vs HLLC vs LLF
}
```

---

## Performance Tuning

### Block Size Optimization

Test different block sizes for your GPU:

```cpp
// Test block sizes: 128, 256, 512
for (int block_size : {128, 256, 512}) {
    auto start = std::chrono::high_resolution_clock::now();
    
    launch_roe_flux(..., block_size);
    cudaDeviceSynchronize();
    
    auto end = std::chrono::high_resolution_clock::now();
    // Measure execution time
}
```

### Memory Access Patterns

If performance is not as expected:

1. **Check coalesced access**: Use `nvprof` or `nsys` to profile
2. **Shared memory**: Consider using shared memory for neighbor data
3. **Register pressure**: Monitor register usage with `--ptxas-options=-v`

Expected block sizes for different GPUs:
- **Ampere (RTX 30xx)**: 256-512 threads/block
- **Turing (RTX 20xx)**: 256 threads/block
- **Pascal (GTX 10xx)**: 128-256 threads/block

---

## Error Handling

Always check CUDA errors:

```cpp
#define CHECK_CUDA_ERROR(call) \
do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA error in %s:%d: %s\n", \
                __FILE__, __LINE__, cudaGetErrorString(err)); \
        exit(EXIT_FAILURE); \
    } \
} while(0)

// Usage:
CHECK_CUDA_ERROR(launch_roe_flux(...));
```

---

## Expected Performance

Based on typical CFD workloads:

| Scheme | Order | Complexity | Expected Speedup (vs CPU) |
|--------|-------|-----------|---------------------------|
| Roe | 1st | High | 25-50x |
| Roe | 2nd | Very High | 20-45x |
| HLLC | 1st | Medium | 30-60x |
| LLF | 1st | Low | 40-70x |
| MUSCL | 2nd | Medium | 20-45x |
| WENO5 | 5th | Very High | 15-35x |

*Note: Actual speedup depends on grid size, GPU model, and memory bandwidth.*

---

## Troubleshooting

### Problem 1: Compilation Errors

**Error**: `undefined reference to roe_flux_kernel`

**Solution**: Ensure wrapper files are included in CMakeLists.txt

### Problem 2: Runtime Errors

**Error**: `illegal memory access`

**Solution**: Check array bounds and neighbor indices
- Verify `face_neighbors` contains valid cell indices
- Check `cell_neighbors` for boundary cells

### Problem 3: Numerical Instability

**Error**: Solution diverges with Roe scheme

**Solution**: 
- Enable entropy fix (already implemented)
- Reduce time step (CFL number)
- Try HLLC or LLF for more robustness

### Problem 4: Performance Lower Than Expected

**Solution**:
- Profile with `nvprof` or Nsight Systems
- Check GPU occupancy
- Optimize block size
- Verify memory coalescing

---

## Next Steps

1. **Immediate**: 
   - ✅ Compile new kernels
   - ✅ Run basic smoke tests
   
2. **Short-term** (this week):
   - Validate against analytical solutions
   - Compare CPU vs GPU results
   - Measure performance speedup
   
3. **Medium-term** (next week):
   - Run full test suite
   - Benchmark different schemes
   - Document best practices
   
4. **Long-term** (next month):
   - Optimize for specific GPU architectures
   - Implement adaptive scheme selection
   - Add 3D support (Priority 3)

---

## References

1. **Roe Scheme**: P.L. Roe, "Approximate Riemann Solvers", JCP, 1981
2. **HLLC Scheme**: E.F. Toro et al., "Restoration of contact surface", Shock Waves, 1994
3. **WENO5**: Jiang & Shu, "Efficient Implementation of WENO", JCP, 1996
4. **MUSCL**: Van Leer, "Towards the ultimate conservative scheme V", JCP, 1979

---

## Contact

For questions or issues with integration:
- Check existing test cases in `Unit_Test_Codes/`
- Review CPU implementations in `src/`
- Consult CUDA documentation in `CUDA_KERNELS/Cuda_Kernels.h`

**Implementation Status**: ✅ COMPLETE - Ready for integration and testing!
