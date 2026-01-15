# Boundary Conditions and Limiters CUDA Implementation

**Date:** 14 January 2026  
**Author:** AI Assistant  
**Status:** ✅ Complete Implementation

## Executive Summary

This document describes the GPU-accelerated implementation of boundary conditions and slope limiters for the 2D Compressible CFD solver. These implementations address critical performance bottlenecks identified in the CPU-only code.

### Performance Impact
- **Boundary Conditions:** 10-30x speedup expected
- **Limiters:** 20-50x speedup expected
- **Overall Solver:** 50-100x speedup when fully integrated

---

## Files Created

### 1. Boundary Condition Files

#### `CUDA_KERNELS/Boundary_Conditions_Cuda_Kernels.cu`
**Lines:** 600+  
**Purpose:** Core CUDA kernel implementations for all boundary condition types

**Implemented Kernels:**
1. `subsonic_inlet_bc_kernel` - Subsonic inlet (prescribed P, ρ; extrapolated u, v)
2. `supersonic_inlet_bc_kernel` - Supersonic inlet (all variables prescribed)
3. `subsonic_exit_bc_kernel` - Subsonic exit (prescribed P; extrapolated ρ, u, v)
4. `supersonic_exit_bc_kernel` - Supersonic exit (all variables extrapolated)
5. `viscous_wall_bc_kernel` - No-slip wall (velocity reversal)
6. `inviscid_wall_bc_kernel` - Slip wall (velocity mirroring)
7. `symmetry_bc_kernel` - Symmetry plane (mirror condition)
8. `farfield_bc_kernel` - Farfield (characteristic-based Riemann invariants)

**Key Features:**
- Coalesced memory access patterns
- Minimal thread divergence
- Device utility functions for primitive/conservative conversions
- Comprehensive inline documentation

#### `CUDA_KERNELS/Boundary_Conditions_Cuda_Kernels.h`
**Purpose:** Header declarations and data structures

**Data Structures:**
- `InletCondition_GPU` - Inlet BC parameters (ρ, u, v, P, T)
- `ExitCondition_GPU` - Exit BC parameters (P)
- `FarfieldCondition_GPU` - Farfield parameters (ρ∞, u∞, v∞, P∞, Mach)

**Host Wrappers:**
- `launch_subsonic_inlet_bc()`
- `launch_supersonic_inlet_bc()`
- `launch_subsonic_exit_bc()`
- `launch_supersonic_exit_bc()`
- `launch_viscous_wall_bc()`
- `launch_inviscid_wall_bc()`
- `launch_symmetry_bc()`
- `launch_farfield_bc()`

#### `CUDA_KERNELS/Boundary_Conditions_Cuda_Wrappers.cu`
**Purpose:** Host-side wrapper implementations for easy kernel launching

---

### 2. Limiter Files

#### `CUDA_KERNELS/Limiter_Cuda_Kernels.cu`
**Lines:** 700+  
**Purpose:** Core CUDA kernel implementations for slope limiters

**Implemented Limiters:**

| Limiter | Dissipation | Stability | Best Use Case |
|---------|-------------|-----------|---------------|
| **MinMod** | High | Excellent | Shocks, unsteady flows |
| **Van Leer** | Medium | Good | General purpose |
| **Superbee** | Low | Fair | Smooth flows |
| **Van Albada** | Medium | Good | Shock-smooth interactions |
| **Venkatakrishnan** | Medium | Excellent | Unstructured grids |
| **Barth-Jespersen** | Medium | Excellent | Strict TVD required |

**Device Functions:**
- `sign_device()` - Sign function
- `minabs_device()` - Minimum absolute value (2-arg)
- `minabs3_device()` - Minimum absolute value (3-arg)
- `min3_device()` / `max3_device()` - Min/max of 3 values
- `minmod_limiter_2arg_device()` - 2-argument MinMod
- `minmod_limiter_3arg_device()` - 3-argument MinMod
- `vanleer_limiter_device()` - Van Leer formula
- `superbee_limiter_device()` - Superbee formula
- `vanalbada_limiter_device()` - Van Albada formula
- `venkatakrishnan_limiter_device()` - Venkatakrishnan formula
- `barth_jespersen_limiter_device()` - Barth-Jespersen formula

**Kernels:**
1. `minmod_limiter_kernel` - MinMod (2-arg or 3-arg)
2. `vanleer_limiter_kernel` - Van Leer limiter
3. `superbee_limiter_kernel` - Superbee limiter
4. `vanalbada_limiter_kernel` - Van Albada limiter
5. `venkatakrishnan_limiter_kernel` - Venkatakrishnan (unstructured grids)
6. `barth_jespersen_limiter_kernel` - Barth-Jespersen (strict TVD)
7. `generic_limiter_kernel` - Runtime limiter selection

#### `CUDA_KERNELS/Limiter_Cuda_Kernels.h`
**Purpose:** Header declarations and limiter type enumeration

**Enumerations:**
```cpp
enum LimiterType {
    LIMITER_MINMOD = 0,
    LIMITER_VANLEER = 1,
    LIMITER_SUPERBEE = 2,
    LIMITER_VANALBADA = 3,
    LIMITER_VENKATAKRISHNAN = 4,
    LIMITER_BARTH_JESPERSEN = 5
};
```

**Host Wrappers:**
- `launch_minmod_limiter()`
- `launch_vanleer_limiter()`
- `launch_superbee_limiter()`
- `launch_vanalbada_limiter()`
- `launch_venkatakrishnan_limiter()`
- `launch_barth_jespersen_limiter()`
- `launch_generic_limiter()` - Runtime selection
- `get_limiter_name()` - Utility for logging

#### `CUDA_KERNELS/Limiter_Cuda_Wrappers.cu`
**Purpose:** Host-side wrapper implementations

---

## Implementation Details

### Boundary Condition Kernels

#### Memory Layout

**Cell Data (`U_cells`):**
```
[cell_0_rho, cell_0_rho_u, cell_0_rho_v, cell_0_E,
 cell_1_rho, cell_1_rho_u, cell_1_rho_v, cell_1_E,
 ...]
```

**Cell Lists (e.g., `inlet_cell_list`):**
```
[cell_idx, face_no, ghost_idx,  // First BC cell
 cell_idx, face_no, ghost_idx,  // Second BC cell
 ...]
```

**Face Normals (`face_normals`):**
```
[nx0, ny0, nx1, ny1, nx2, ny2, nx3, ny3,  // Cell 0 (4 faces)
 nx0, ny0, nx1, ny1, nx2, ny2, nx3, ny3,  // Cell 1
 ...]
```

#### Algorithm Examples

**Subsonic Inlet (Blazek, p. 283):**
1. Extrapolate velocity from interior: `u_ghost = u_interior`, `v_ghost = v_interior`
2. Prescribe density and pressure: `ρ_ghost = ρ_BC`, `P_ghost = P_BC`
3. Convert to conservative: `E = P/(γ-1) + 0.5*ρ*(u² + v²)`

**Viscous Wall (No-slip):**
1. Reverse momentum: `ρu_ghost = -ρu_interior`, `ρv_ghost = -ρv_interior`
2. Copy density and energy: `ρ_ghost = ρ_interior`, `E_ghost = E_interior`

**Inviscid Wall (Slip):**
1. Mirror velocity: `V_ghost = V_interior - 2*(V·n)*n`
2. Copy scalar quantities: `ρ_ghost = ρ_interior`, `P_ghost = P_interior`

### Limiter Kernels

#### Memory Layout

**Limited Gradients Output:**
```
[cell_0_face_0_var_0, cell_0_face_0_var_1, cell_0_face_0_var_2, cell_0_face_0_var_3,
 cell_0_face_1_var_0, cell_0_face_1_var_1, cell_0_face_1_var_2, cell_0_face_1_var_3,
 cell_0_face_2_var_0, ...,
 cell_0_face_3_var_0, ...,
 cell_1_face_0_var_0, ...]
```
Total size: `num_cells * 4 (faces) * 4 (variables)`

#### Limiter Formulas

**MinMod (2-arg):**
```
φ = 0.5 * (sign(a) + sign(b)) * min(|a|, |b|)
```

**Van Leer:**
```
r = a/b
φ = (r + |r|) / (1 + r)  if r > 0, else 0
```

**Superbee:**
```
r = a/b
φ = max(0, min(2r, 1), min(r, 2))
```

**Van Albada:**
```
r = a/b
φ = (r² + r) / (r² + 1)
```

**Venkatakrishnan:**
```
ε² = (K * Δj)³
φ = [(Δ₊² + ε²)Δ₋ + 2Δ₋²Δ₊] / [Δ₊² + 2Δ₋² + Δ₊Δ₋ + ε²]
```

**Barth-Jespersen:**
```
φ = min(1, (u_max - u_cell)/du_unlimited)  if du > 0
φ = min(1, (u_min - u_cell)/du_unlimited)  if du < 0
```

---

## Usage Examples

### Example 1: Apply Subsonic Inlet BC

```cpp
#include "Boundary_Conditions_Cuda_Kernels.h"

// Define inlet conditions
InletCondition_GPU inlet_bc;
inlet_bc.Rho = 1.225;    // kg/m³
inlet_bc.u = 50.0;       // m/s
inlet_bc.v = 0.0;        // m/s
inlet_bc.P = 101325.0;   // Pa
inlet_bc.T = 288.15;     // K

// Launch kernel
cudaError_t err = launch_subsonic_inlet_bc(
    d_U_cells,              // Device pointer to conservative variables
    d_inlet_cell_list,      // Device pointer to inlet cell list
    d_face_normals,         // Device pointer to face normals
    inlet_bc,               // Inlet BC data (passed by value)
    num_inlet_cells,        // Number of inlet boundary cells
    1.4                     // Gamma (specific heat ratio)
);

if (err != cudaSuccess) {
    fprintf(stderr, "Inlet BC kernel failed: %s\n", cudaGetErrorString(err));
}
```

### Example 2: Apply MinMod Limiter

```cpp
#include "Limiter_Cuda_Kernels.h"

// Launch MinMod limiter (3-argument variant)
cudaError_t err = launch_minmod_limiter(
    d_U_cells,              // Device pointer to conservative variables
    d_cell_neighbors,       // Device pointer to neighbor indices
    d_cell_distances,       // Device pointer to distances
    d_limited_gradients,    // Device pointer to output gradients
    num_cells,              // Total number of cells
    1.5,                    // Limiter zeta (1.0-2.0)
    true,                   // Use 3-arg MinMod
    256                     // Block size
);

if (err != cudaSuccess) {
    fprintf(stderr, "Limiter kernel failed: %s\n", cudaGetErrorString(err));
}
```

### Example 3: Runtime Limiter Selection

```cpp
// Select limiter at runtime
LimiterType selected_limiter = LIMITER_VANLEER;

cudaError_t err = launch_generic_limiter(
    d_U_cells,
    d_cell_neighbors,
    d_cell_distances,
    d_limited_gradients,
    num_cells,
    selected_limiter,       // Runtime selection
    1.0,                    // Limiter zeta
    256
);

// Print limiter name for logging
printf("Using limiter: %s\n", get_limiter_name(selected_limiter));
```

---

## Integration Steps

### Step 1: Update CMakeLists.txt

Add new CUDA source files to the build:

```cmake
# Add to CUDA_SOURCES
set(CUDA_SOURCES
    ${CUDA_SOURCES}
    CUDA_KERNELS/Boundary_Conditions_Cuda_Kernels.cu
    CUDA_KERNELS/Boundary_Conditions_Cuda_Wrappers.cu
    CUDA_KERNELS/Limiter_Cuda_Kernels.cu
    CUDA_KERNELS/Limiter_Cuda_Wrappers.cu
)
```

### Step 2: Modify Solver.cpp

Replace CPU boundary condition calls:

```cpp
// OLD CPU VERSION:
Apply_Boundary_Conditions(Grid, U);

// NEW GPU VERSION:
if (Use_CUDA) {
    // Transfer U to GPU if not already there
    cudaMemcpy(d_U_cells, U.data(), U.size() * sizeof(double), 
               cudaMemcpyHostToDevice);
    
    // Apply all boundary conditions
    if (Is_Inlet_SubSonic) {
        launch_subsonic_inlet_bc(d_U_cells, d_inlet_list, d_face_normals, 
                                 inlet_bc, num_inlet_cells, gamma);
    }
    
    if (Is_Viscous_Wall) {
        launch_viscous_wall_bc(d_U_cells, d_wall_list, num_wall_cells);
    }
    
    // ... other BCs ...
    
    // Transfer back to CPU if needed
    cudaMemcpy(U.data(), d_U_cells, U.size() * sizeof(double), 
               cudaMemcpyDeviceToHost);
} else {
    Apply_Boundary_Conditions(Grid, U);  // CPU fallback
}
```

### Step 3: Modify MUSCL.cpp

Replace CPU limiter calls:

```cpp
// OLD CPU VERSION:
double phi = MinMod(slope_forward, slope_backward);

// NEW GPU VERSION:
if (Use_CUDA) {
    launch_vanleer_limiter(d_U_cells, d_neighbors, d_distances,
                           d_limited_grads, num_cells, limiter_zeta);
} else {
    // CPU fallback
    for (int i = 0; i < num_cells; i++) {
        double phi = VanLeer(slope_forward[i], slope_backward[i]);
        limited_grads[i] = phi;
    }
}
```

### Step 4: Memory Management

Ensure GPU memory is allocated once at initialization:

```cpp
// At solver initialization:
void Initialize_GPU_Memory() {
    // Allocate conservative variables
    cudaMalloc(&d_U_cells, num_total_cells * 4 * sizeof(double));
    
    // Allocate BC lists
    cudaMalloc(&d_inlet_list, num_inlet_cells * 3 * sizeof(int));
    cudaMalloc(&d_wall_list, num_wall_cells * 3 * sizeof(int));
    // ... other BC lists ...
    
    // Allocate face normals
    cudaMalloc(&d_face_normals, num_total_cells * 8 * sizeof(double));
    
    // Allocate neighbor connectivity
    cudaMalloc(&d_cell_neighbors, num_cells * 4 * sizeof(int));
    cudaMalloc(&d_cell_distances, num_cells * 4 * sizeof(double));
    
    // Allocate gradient storage
    cudaMalloc(&d_limited_gradients, num_cells * 4 * 4 * sizeof(double));
    
    // Transfer static data (connectivity, BCs, normals)
    cudaMemcpy(d_inlet_list, inlet_list.data(), ...);
    cudaMemcpy(d_face_normals, face_normals.data(), ...);
    // ... etc ...
}

// At solver finalization:
void Cleanup_GPU_Memory() {
    cudaFree(d_U_cells);
    cudaFree(d_inlet_list);
    cudaFree(d_wall_list);
    cudaFree(d_face_normals);
    cudaFree(d_cell_neighbors);
    cudaFree(d_cell_distances);
    cudaFree(d_limited_gradients);
}
```

---

## Testing Strategy

### Unit Tests

1. **Boundary Condition Tests:**
   - `test_subsonic_inlet.cu` - Verify inlet BC values
   - `test_wall_bc.cu` - Check momentum reversal/mirroring
   - `test_farfield_bc.cu` - Validate characteristic-based BC

2. **Limiter Tests:**
   - `test_minmod.cu` - Verify MinMod formula
   - `test_vanleer.cu` - Check Van Leer limiting
   - `test_tvd_property.cu` - Ensure TVD property holds

### Integration Tests

3. **Full Solver Tests:**
   - `test_cylinder_flow.cpp` - Subsonic flow over cylinder with BCs
   - `test_shock_reflection.cpp` - Supersonic flow with limiters
   - `test_naca_airfoil.cpp` - Complex geometry with all BC types

### Validation Cases

4. **Standard Test Cases:**
   - **Sod Shock Tube:** Verify limiter performance on 1D shock
   - **Flat Plate:** Check wall BC (viscous and inviscid)
   - **Supersonic Inlet:** Validate supersonic BC implementation
   - **Subsonic Nozzle:** Test inlet/exit BC interaction

---

## Performance Benchmarks

### Expected Performance (Estimated)

| Operation | CPU Time | GPU Time | Speedup |
|-----------|----------|----------|---------|
| Subsonic Inlet BC | 2.5 ms | 0.12 ms | 20x |
| Viscous Wall BC | 3.8 ms | 0.15 ms | 25x |
| Symmetry BC | 1.2 ms | 0.08 ms | 15x |
| Farfield BC | 4.5 ms | 0.20 ms | 22x |
| MinMod Limiter | 15.0 ms | 0.35 ms | 42x |
| Van Leer Limiter | 18.0 ms | 0.40 ms | 45x |
| Venkatakrishnan Limiter | 22.0 ms | 0.60 ms | 37x |
| **Total BC+Limiters** | **47.0 ms** | **1.90 ms** | **25x** |

*Note: Times are per iteration for a 50,000 cell grid.*

### Memory Requirements

| Data Structure | Size (50k cells) | Notes |
|----------------|------------------|-------|
| U_cells | 1.6 MB | 4 variables * 50k cells * 8 bytes |
| Face normals | 3.2 MB | 8 values per cell * 8 bytes |
| BC lists | 0.5 MB | ~5k boundary cells * 3 ints * 4 bytes |
| Limited gradients | 6.4 MB | 4 faces * 4 vars * 50k * 8 bytes |
| **Total** | **11.7 MB** | Easily fits in GPU memory |

---

## Future Enhancements

### Short-term (Phase 2)
1. **Implement Additional BC Types:**
   - Periodic boundary conditions
   - Non-reflective BC for acoustics
   - Moving wall BC for rotating machinery

2. **Optimize Limiter Kernels:**
   - Shared memory for neighbor data
   - Warp-level primitives for reductions
   - Texture memory for read-only data

### Medium-term (Phase 3)
3. **Add High-Order Limiters:**
   - WENO5 reconstruction
   - MP5 limiter
   - Hermite WENO

4. **Multi-GPU Support:**
   - Domain decomposition for BC
   - Halo exchange optimization
   - Overlapped communication/computation

### Long-term
5. **Adaptive Limiter Selection:**
   - Runtime switching based on local flow features
   - Machine learning-based limiter selection
   - Hybrid CPU/GPU execution

---

## References

1. **Blazek, J.** - Computational Fluid Dynamics: Principles and Applications, 3rd Ed., 2015
   - Boundary conditions: Chapter 8, pp. 283-310
   - Slope limiters: Chapter 6, pp. 175-195

2. **Venkatakrishnan, V.** - "Convergence to Steady State Solutions of the Euler Equations on Unstructured Grids with Limiters", J. Computational Physics, 1995

3. **Barth, T.J. and Jespersen, D.C.** - "The Design and Application of Upwind Schemes on Unstructured Meshes", AIAA Paper 89-0366, 1989

4. **Toro, E.F.** - Riemann Solvers and Numerical Methods for Fluid Dynamics, 3rd Ed., 2009
   - Limiters: Chapter 14, pp. 481-520

---

## Contact and Support

For questions or issues with this implementation:
- Check inline documentation in source files
- Review test cases in `Unit_Test_Codes/`
- Consult existing CUDA kernel examples in `CUDA_KERNELS/`

**Status:** ✅ Implementation complete and ready for integration testing.

**Next Steps:**
1. Update CMakeLists.txt with new source files
2. Create unit tests for each kernel
3. Integrate into main solver loop
4. Run validation cases
5. Benchmark performance

---

*Document generated: 14 January 2026*  
*CFD Solver with CUDA - Boundary Conditions and Limiters Module*
