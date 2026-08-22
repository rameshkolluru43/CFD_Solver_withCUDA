# Configuration Reference

The solver is driven by **JSON** configuration files. The main entry is typically a top-level run config that references a test-case JSON.

## Top-Level Structure

```json
{
  "TestCase": { ... },
  "Simulation": { ... },
  "Solver": { ... },
  "LimiterCoefficients": { ... }
}
```

## TestCase

| Key | Type | Description |
|-----|------|-------------|
| `Test_Case` | int | Test case identifier. |
| `Test_Case_Name` | string | Name (e.g. `"Half_Cylinder"`). |
| `Test_Case_Json` | string | Path to test case JSON (e.g. `"../json_Files/Half_Cylinder_481_161_M6_viscous.json"`). |

The test-case JSON defines grid path, boundary conditions, Re / Tw, and flow parameters.

## Simulation

| Key | Type | Description |
|-----|------|-------------|
| `Initialize_Type` | int | `0` freestream/IC from case; `1` restart / continue from prior solution (used by corner-fix smoke). |
| `Is_Implicit_Method` | bool | Implicit time stepping (stub in main loop — leave `false`). |
| `Total_Iterations` | int | Maximum iterations. |
| `CFL` | double | CFL number (viscous M=6 continues often use `0.02`). |
| `Terminating_Time` | double | Stop time for time-dependent runs. |
| `Is_Time_Dependent` | bool | Time-accurate vs steady. |

## Solver

| Key | Type | Description |
|-----|------|-------------|
| `Solver_Type` | int | Solver type. |
| `Solver_Name` | string | e.g. `"Euler"`, `"Navier_Stokes"`. |
| `Is_Conservative` | bool | Conservative formulation. |
| `Is_Viscous` | bool | Enable Navier–Stokes / viscous wall path (`Is_Viscous_Wall`). |
| `Limiter_Case` | int | Limiter choice. |
| `Area_Weighted_Average` | int | Gradient/weighting option. |
| `Flux_Type` | int | **Output directory label only** — physics uses `Dissipation_Type`. |
| `NUM_FLUX_COMPONENTS` | int | Usually 4 (2D Euler/NS). |
| `Is_Second_Order` | bool | MUSCL second-order spatial (ignored when `Is_WENO`). |
| `Time_Accurate` | bool | TVD-RK vs local/explicit. |
| `Local_Time_Stepping` | bool | Local time stepping. |
| `Non_Dimensional_Form` | bool | Non-dimensional form. |
| `Is_WENO` | bool | Host WENO reconstruction (quad meshes). |
| `Is_Char` | bool | Characteristic WENO (when supported). |
| `Dissipation_Type` | int | See table below. |
| `Is_MOVERS_1` | bool | MOVERS variant. |
| `Enable_Entropy_Fix` | bool | Roe / scheme entropy fix. |
| **AMR** | | |
| `Enable_AMR` | bool | Enable gradient-based AMR (inviscid loop). |
| `AMR_Period` | int | Apply AMR every N iterations. |
| `AMR_Start_Iteration` | int | First iteration for AMR. |
| `AMR_Gradient_Threshold` | double | Refine if indicator > this. |
| `AMR_Max_Fraction` | double | Max fraction of cells to tag (0 or ≥1 = no cap). |
| `AMR_Coarsen_Threshold` | double | Optional coarsening threshold. |

### Dissipation_Type map

| Value | Scheme | Notes |
|------:|--------|-------|
| 1 | LLF | CUDA 1O main-loop |
| 2 | MOVERS | CUDA 1O main-loop |
| 3 | Roe | Host (CUDA dispatcher falls back) |
| 4 | RICCA | CUDA 1O; **WENO face dissipation** on host |
| 5 | MOVERS_NWSC | CUDA 1O main-loop |
| 6 | RICCA_LLF hybrid | WENO pressure-sensor blend toward LLF |

## Test-case JSON (flow / wall)

Keys live under the file pointed to by `Test_Case_Json` (e.g. `Half_Cylinder_481_161_M6_viscous.json`):

| Area | Typical keys |
|------|----------------|
| Mesh | `mesh_parameters.Nx`, `Ny`, `Grid_Size`, grid path via test-case loader |
| Freestream / inlet | `Flow_Conditions.inlet_conditions.*` (Mach, ρ, p, u, v) |
| Wall | `wall_type` (`no_slip`), `wall_temperature` (Tw), `wall_velocity` |
| Viscous | `Reynolds_number` (Re), `Prandtl_number` |

Example viscous half-cylinder: Re = `1e5`, Tw = `2.0`, M∞ = `6`, isothermal no-slip.

## LimiterCoefficients

| Key | Type | Description |
|-----|------|-------------|
| `Limiter_Zeta` | double | Limiter parameter (e.g. 0.5). |
| `Limiter_Zeta1` | double | Second limiter parameter. |
| `Venkat_K` | double | Venkatakrishnan parameter (when used). |

## Example: half-cylinder viscous smoke (corner fix)

See `input/json_Files/Run_HalfCylinder_M6_corner_fix_smoke500.json`:

```json
{
  "TestCase": {
    "Test_Case": 1,
    "Test_Case_Name": "Half_Cylinder",
    "Test_Case_Json": "../json_Files/Half_Cylinder_481_161_M6_viscous.json"
  },
  "Simulation": {
    "Initialize_Type": 1,
    "Is_Implicit_Method": false,
    "Total_Iterations": 500,
    "CFL": 0.02,
    "Terminating_Time": 1000.0,
    "Is_Time_Dependent": false
  },
  "Solver": {
    "Solver_Type": 0,
    "Solver_Name": "Navier_Stokes",
    "Is_Conservative": true,
    "Is_Viscous": true,
    "Is_Second_Order": false,
    "Local_Time_Stepping": true,
    "Non_Dimensional_Form": true,
    "Is_WENO": true,
    "Dissipation_Type": 4,
    "Enable_Entropy_Fix": true
  }
}
```

Full validation recipe: [HALF_CYLINDER_VALIDATION.md](HALF_CYLINDER_VALIDATION.md).

## Example: AMR Enabled

See `input/json_Files/Test_Config_AMR.json` for gradient AMR tagging on the inviscid path.

## How Config Is Loaded

- Command-line argument is the run JSON path (e.g. `input/json_Files/Run_HalfCylinder_M6_corner_fix_smoke500.json`).
- **Configuration_Read** parses `TestCase`, `Simulation`, `Solver`, and `LimiterCoefficients`.
- Test-case JSON supplies grid and BC / Re / Tw.
- AMR keys are optional (defaults: `Enable_AMR=false`, etc.).
