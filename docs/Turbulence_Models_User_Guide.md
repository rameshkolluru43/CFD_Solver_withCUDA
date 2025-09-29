# Turbulence Models Implementation Guide

## Overview

This implementation adds K-epsilon and K-omega turbulence models to the 2D compressible CFD solver. The implementation includes:

- **K-epsilon turbulence model** (Standard Launder-Sharma)
- **K-omega turbulence model** (Standard Wilcox)  
- **SST K-omega turbulence model** (Menter)
- **Wall functions** for high Reynolds number flows
- **Low-Re treatment** for near-wall resolution
- **Validation test cases** for model verification

## File Structure

### Header Files
- `include/Turbulence_Models.h` - Main turbulence model declarations
- Contains constants, structures, and function prototypes

### Implementation Files
- `src/Turbulence_Models.cpp` - Main turbulence model interface
- `src/K_Epsilon_Model.cpp` - K-epsilon model implementation
- `src/K_Omega_Model.cpp` - K-omega and SST model implementation
- `src/Turbulence_Integration.cpp` - Integration with main solver

### Test Cases
- `Test_Cases/Turbulent_Flat_Plate_Validation.cpp` - Validation test case
- Implements turbulent boundary layer over flat plate

### Configuration
- `json_Files/Turbulent_Flow_2D.json` - Example configuration file
- Contains all turbulence model parameters and settings

## Usage Instructions

### 1. Basic Setup

```cpp
#include "Turbulence_Models.h"

// Initialize turbulence model
TurbulenceModel model = TurbulenceModel::K_EPSILON;
Initialize_Turbulence_Model(model);

// Set characteristic length scale
characteristic_length = 1.0;  // meters

// Enable wall functions
use_wall_functions = true;
y_plus_target = 30.0;
```

### 2. Configuration File Setup

Create a JSON configuration file (see `Turbulent_Flow_2D.json`):

```json
{
  "turbulence_model": {
    "model_type": "k_epsilon",
    "use_wall_functions": true,
    "y_plus_target": 30.0,
    "turbulence_intensity": 0.05,
    "turbulent_length_scale_ratio": 0.1
  }
}
```

### 3. Integration with Main Solver

Replace the standard viscous solver call:

```cpp
// Instead of:
// Viscous_Solver(error_file, solution_file);

// Use:
Viscous_Solver_With_Turbulence(error_file, solution_file);
```

### 4. Boundary Conditions

Set appropriate turbulence boundary conditions:

```cpp
// Inlet conditions
Apply_Turbulence_Inlet_BC(cell_index);

// Wall conditions  
Apply_Turbulence_Wall_BC(cell_index);

// Outlet conditions
Apply_Turbulence_Outlet_BC(cell_index);
```

## Turbulence Models

### K-epsilon Model

The standard K-epsilon model equations:

**k-equation:**
```
∂k/∂t + ∇·(ρuk) = ∇·[(μ + μt/σk)∇k] + Pk - ρε
```

**ε-equation:**
```
∂ε/∂t + ∇·(ρuε) = ∇·[(μ + μt/σε)∇ε] + C1(ε/k)Pk - C2ρ(ε²/k)
```

**Turbulent viscosity:**
```
μt = ρCμk²/ε
```

**Model constants:**
- Cμ = 0.09
- C1 = 1.44  
- C2 = 1.92
- σk = 1.0
- σε = 1.3

### K-omega Model

The standard K-omega model equations:

**k-equation:**
```
∂k/∂t + ∇·(ρuk) = ∇·[(μ + σ*μt)∇k] + Pk - β*ρkω
```

**ω-equation:**
```
∂ω/∂t + ∇·(ρuω) = ∇·[(μ + σμt)∇ω] + α(ω/k)Pk - βρω²
```

**Turbulent viscosity:**
```
μt = ρk/ω
```

**Model constants:**
- β* = 0.09
- β = 0.075
- σ = 0.5
- σ* = 0.5  
- α = 5/9

### SST K-omega Model

The SST model blends K-omega and K-epsilon formulations:

**Blending function F1:**
```
F1 = tanh(arg1⁴)
```

**Turbulent viscosity:**
```
μt = ρa₁k/max(a₁ω, ΩF₂)
```

Where Ω is the vorticity magnitude and F₂ is the second blending function.

## Wall Treatment

### Wall Functions (y⁺ > 30)

For high Reynolds number flows, wall functions are used:

**Velocity profile:**
```
u⁺ = (1/κ)ln(y⁺) + B
```

**Turbulence variables at wall:**
```
k = u_τ²/√Cμ
ε = u_τ³/(κy)    (K-epsilon)
ω = u_τ/(√β*κy) (K-omega)
```

### Low-Re Treatment (y⁺ < 1)

For low Reynolds number flows:

```
k_wall = 0
ε_wall = 2μk/(ρy²)     (K-epsilon)
ω_wall = 60μ/(ρβy²)    (K-omega)
```

## Boundary Conditions

### Inlet Boundary

Set turbulence intensity and length scale:

```cpp
double Tu = 0.05;  // 5% turbulence intensity
double L = 0.1 * characteristic_length;  // Length scale

double k_inlet = 1.5 * (Tu * U_inlet)²
double ε_inlet = Cμ * k^1.5 / L        (K-epsilon)  
double ω_inlet = k / (√β* * L)          (K-omega)
```

### Wall Boundary

- **No-slip condition** for velocity
- **Wall functions** or **low-Re treatment** for turbulence
- **Adiabatic** or **isothermal** wall temperature

### Outlet Boundary

- **Zero gradient** for all turbulence variables
- **Extrapolation** from upstream cells

### Symmetry Boundary

- **Zero gradient** for k, ε, ω
- **Zero normal velocity** component

## Validation and Testing

### Flat Plate Validation

Run the flat plate validation test:

```cpp
// Test K-epsilon model
bool ke_passed = Run_Turbulent_Flat_Plate_Validation(TurbulenceModel::K_EPSILON);

// Test K-omega model  
bool ko_passed = Run_Turbulent_Flat_Plate_Validation(TurbulenceModel::K_OMEGA_SST);
```

**Expected results:**
- Skin friction coefficient: Cf = 0.00332 ± 10%
- Boundary layer thickness: δ/L = 0.382/Re_L^(1/5)
- Displacement thickness: δ*/L = 0.0463/Re_L^(1/5)

### Other Test Cases

Additional validation cases can be implemented:

1. **Channel Flow** - Fully developed turbulent channel flow
2. **Backward Facing Step** - Separated flow reattachment
3. **Flow Over Cylinder** - Bluff body turbulent wake

## Output and Post-Processing

### Turbulence Variables Output

The solver outputs the following turbulence variables:

```
Cell_Index X Y Z k epsilon omega mu_t y_plus
```

### Visualization

For VTK output, turbulence variables are included:

- Turbulent kinetic energy (k)
- Turbulent dissipation rate (ε) 
- Specific dissipation rate (ω)
- Turbulent viscosity (μt)
- Y-plus values (y⁺)

### Convergence Monitoring

Monitor convergence of turbulence equations:

```
Iter    dt        Rho_Error   k_Error     eps_Error   Wall_Time
1000    1.23e-5   2.45e-4     1.23e-5     3.45e-6     0.234
```

## Performance Considerations

### Grid Requirements

**K-epsilon model:**
- Can use wall functions with y⁺ = 30-300
- Coarser near-wall grids acceptable

**K-omega/SST model:**
- Better near-wall performance
- Requires y⁺ < 1 for best results
- More sensitive to grid quality

### Computational Cost

- Adds ~40-50% to computational cost
- Two additional transport equations
- Additional source terms and diffusion

### Memory Requirements

- Additional storage for k, ε/ω variables
- Production terms and gradients
- Wall distance calculations

## Troubleshooting

### Common Issues

1. **Negative turbulence variables**
   - Check for proper limiters
   - Verify boundary conditions
   - Reduce time step if needed

2. **Excessive turbulent viscosity**
   - Apply realizability constraints
   - Check production limiters
   - Verify model constants

3. **Poor convergence**
   - Reduce time step
   - Use implicit treatment
   - Check grid quality

4. **Unphysical results**
   - Verify wall distance calculation
   - Check boundary condition implementation
   - Validate against test cases

### Debugging Tips

1. **Enable detailed output:**
   ```cpp
   Calculate_Turbulence_Statistics();
   Write_Turbulence_Variables("debug_turb.dat");
   ```

2. **Check y⁺ distribution:**
   ```cpp
   // Values should be y⁺ > 30 for wall functions
   // or y⁺ < 1 for low-Re treatment
   ```

3. **Monitor production terms:**
   ```cpp
   // Production should be positive and reasonable
   // Pk ~ O(k*ε/k) for k-epsilon
   ```

## References

1. Launder, B.E. and Sharma, B.I. (1974) - K-epsilon model
2. Wilcox, D.C. (2006) - K-omega model  
3. Menter, F.R. (1994) - SST K-omega model
4. Pope, S.B. (2000) - Turbulent Flows
5. Blazek, J. (2015) - Computational Fluid Dynamics

## License

This implementation is part of the CFD_Solver_withCUDA project and follows the same license terms.

## Contributing

For contributions or bug reports, please follow the project guidelines and ensure all validation tests pass before submitting changes.