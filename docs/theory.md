# CFD Theory and Mathematical Background

## Overview

This document provides the mathematical foundation and theoretical background for the CFD Solver with CUDA. The solver implements numerical methods for solving the Navier-Stokes equations using finite volume, finite element, or finite difference approaches.

## Governing Equations

### Navier-Stokes Equations

The incompressible Navier-Stokes equations in conservation form are:

**Continuity Equation:**
```
∂ρ/∂t + ∇·(ρu) = 0
```

For incompressible flow (ρ = constant):
```
∇·u = 0
```

**Momentum Equations:**
```
∂(ρu)/∂t + ∇·(ρuu) = -∇p + μ∇²u + ρg
```

Where:
- **u** = velocity vector (u, v, w)
- **p** = pressure
- **ρ** = density
- **μ** = dynamic viscosity
- **g** = gravitational acceleration

### Dimensionless Form

The equations can be made dimensionless using characteristic scales:
- Length scale: L
- Velocity scale: U
- Time scale: L/U
- Pressure scale: ρU²

This leads to the dimensionless parameters:
- **Reynolds number**: Re = ρUL/μ
- **Mach number**: Ma = U/c (where c is speed of sound)

## Numerical Methods

### Finite Volume Method

The finite volume method integrates the governing equations over control volumes:

```
∫∫∫_V (∂φ/∂t) dV + ∫∫_S (φu·n) dS = ∫∫_S (Γ∇φ·n) dS + ∫∫∫_V S_φ dV
```

Where:
- φ = conserved variable
- V = control volume
- S = surface area
- Γ = diffusion coefficient
- S_φ = source term
- n = outward normal vector

### Spatial Discretization

#### Central Differencing
For second-order accuracy:
```
(∂φ/∂x)_i = (φ_{i+1} - φ_{i-1})/(2Δx)
```

#### Upwind Differencing
For convection-dominated flows:
```
(∂φ/∂x)_i = (φ_i - φ_{i-1})/Δx  (for u > 0)
(∂φ/∂x)_i = (φ_{i+1} - φ_i)/Δx  (for u < 0)
```

### Time Integration

#### Explicit Euler
```
φ^{n+1} = φ^n + Δt · R(φ^n)
```

#### Runge-Kutta 4th Order
```
k₁ = Δt · R(φⁿ)
k₂ = Δt · R(φⁿ + k₁/2)
k₃ = Δt · R(φⁿ + k₂/2)
k₄ = Δt · R(φⁿ + k₃)
φ^{n+1} = φⁿ + (k₁ + 2k₂ + 2k₃ + k₄)/6
```

## Pressure-Velocity Coupling

### SIMPLE Algorithm

The Semi-Implicit Method for Pressure-Linked Equations (SIMPLE):

1. **Guess pressure field** p*
2. **Solve momentum equations** to get u*, v*, w*
3. **Solve pressure correction equation**:
   ```
   ∇²p' = ∇·u*/Δt
   ```
4. **Update pressure**: p = p* + p'
5. **Update velocities**: u = u* - Δt(∇p')
6. **Check convergence**, if not converged, go to step 1

### Pressure Poisson Equation

For incompressible flow, the pressure Poisson equation is:
```
∇²p = -ρ(∇·u)/Δt + ρ∇·(u·∇u) - ρ∇·g
```

## Boundary Conditions

### Dirichlet Boundary Conditions
Specify the value of the variable at the boundary:
```
φ_boundary = φ_specified
```

### Neumann Boundary Conditions
Specify the gradient at the boundary:
```
∂φ/∂n|_boundary = q_specified
```

### Wall Boundary Conditions
- **No-slip**: u = v = w = 0
- **Free-slip**: u·n = 0, ∂u_t/∂n = 0

### Inlet/Outlet Conditions
- **Inlet**: Specify velocity profile
- **Outlet**: ∂u/∂n = 0, p = p_outlet

## Turbulence Modeling

### Reynolds-Averaged Navier-Stokes (RANS)

The time-averaged momentum equation:
```
∂(ρŪ)/∂t + ∇·(ρŪŪ) = -∇P̄ + μ∇²Ū - ∇·(ρu'u')
```

Where the Reynolds stress tensor is:
```
τᵢⱼ = -ρu'ᵢu'ⱼ
```

### k-ε Model

Transport equation for turbulent kinetic energy:
```
∂(ρk)/∂t + ∇·(ρUk) = ∇·[(μ + μₜ/σₖ)∇k] + Pₖ - ρε
```

Transport equation for dissipation rate:
```
∂(ρε)/∂t + ∇·(ρUε) = ∇·[(μ + μₜ/σₑ)∇ε] + C₁(ε/k)Pₖ - C₂ρ(ε²/k)
```

Turbulent viscosity:
```
μₜ = ρCμk²/ε
```

## Stability and Convergence

### CFL Condition
For explicit time stepping:
```
CFL = |u|Δt/Δx + |v|Δt/Δy + |w|Δt/Δz ≤ CFL_max
```

Typically CFL_max ≤ 1 for stability.

### Diffusion Number
```
D = νΔt/Δx² ≤ 0.5
```

### Convergence Criteria
Monitor residuals:
```
R = ||φⁿ⁺¹ - φⁿ||₂/||φⁿ||₂ < tolerance
```

## CUDA Implementation Considerations

### Memory Coalescing
- Arrange data structures to enable coalesced memory access
- Use structure of arrays (SoA) instead of array of structures (AoS)

### Thread Block Organization
- Optimal thread block size: 128-512 threads
- Consider shared memory usage
- Minimize warp divergence

### Numerical Precision
- Use double precision for better accuracy
- Consider mixed precision for performance

## References

1. Ferziger, J.H. and Perić, M., "Computational Methods for Fluid Dynamics"
2. Versteeg, H.K. and Malalasekera, W., "An Introduction to Computational Fluid Dynamics"
3. NVIDIA CUDA Programming Guide
4. Blazek, J., "Computational Fluid Dynamics: Principles and Applications"