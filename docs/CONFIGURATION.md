# Configuration Reference

The solver is driven by **JSON** configuration files. The main entry is typically a top-level config that references a test case JSON.

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
| `Test_Case_Json` | string | Path to test case JSON (e.g. `"../json_Files/Half_Cylinder.json"`). |

The test case JSON usually defines grid path, boundary conditions, and flow parameters.

## Simulation

| Key | Type | Description |
|-----|------|-------------|
| `Initialize_Type` | int | Initialization method. |
| `Is_Implicit_Method` | bool | Use implicit time stepping. |
| `Total_Iterations` | int | Maximum iterations. |
| `CFL` | double | CFL number. |
| `Terminating_Time` | double | Stop time for time-dependent runs. |
| `Is_Time_Dependent` | bool | Time-accurate vs steady. |

## Solver

| Key | Type | Description |
|-----|------|-------------|
| `Solver_Type` | int | Solver type. |
| `Solver_Name` | string | e.g. `"Euler"`. |
| `Is_Conservative` | bool | Conservative formulation. |
| `Is_Viscous` | bool | Viscous terms on/off. |
| `Limiter_Case` | int | Limiter choice. |
| `Area_Weighted_Average` | int | Gradient/weighting option. |
| `Flux_Type` | int | Flux scheme selector. |
| `NUM_FLUX_COMPONENTS` | int | Usually 4 (2D Euler). |
| `Is_Second_Order` | bool | Second-order spatial. |
| `Time_Accurate` | bool | Time-accurate integration. |
| `Local_Time_Stepping` | bool | Local time stepping. |
| `Non_Dimensional_Form` | bool | Non-dimensional form. |
| `Is_WENO` | bool | Use WENO reconstruction (quad meshes). |
| `Dissipation_Type` | int | 1=LLF, 2=MOVERS, 3=Roe, 4=RICCA, 5=MOVERS_NWSC. |
| `Is_MOVERS_1` | bool | MOVERS variant. |
| `Enable_Entropy_Fix` | bool | Roe entropy fix. |
| **AMR** | | |
| `Enable_AMR` | bool | Enable gradient-based AMR tagging. |
| `AMR_Period` | int | Apply AMR every N iterations. |
| `AMR_Gradient_Threshold` | double | Refine if indicator > this. |
| `AMR_Max_Fraction` | double | Max fraction of cells to tag (0 or ≥1 = no cap). |

## LimiterCoefficients

| Key | Type | Description |
|-----|------|-------------|
| `Limiter_Zeta` | double | Limiter parameter (e.g. 0.5). |
| `Limiter_Zeta1` | double | Second limiter parameter. |

## Example: AMR Enabled

See `json_Files/Test_Config_AMR.json`:

```json
{
  "TestCase": {
    "Test_Case": 1,
    "Test_Case_Name": "Half_Cylinder",
    "Test_Case_Json": "../json_Files/Half_Cylinder.json"
  },
  "Simulation": {
    "Initialize_Type": 0,
    "Is_Implicit_Method": false,
    "Total_Iterations": 500,
    "CFL": 0.3,
    "Terminating_Time": 1000.0,
    "Is_Time_Dependent": false
  },
  "Solver": {
    "Solver_Type": 0,
    "Solver_Name": "Euler",
    "Is_Conservative": true,
    "Is_Viscous": false,
    "Is_Second_Order": true,
    "Is_WENO": false,
    "Dissipation_Type": 1,
    "Enable_Entropy_Fix": false,
    "Enable_AMR": true,
    "AMR_Period": 100,
    "AMR_Gradient_Threshold": 0.05,
    "AMR_Max_Fraction": 0.25
  },
  "LimiterCoefficients": {
    "Limiter_Zeta": 0.5,
    "Limiter_Zeta1": 0.5
  }
}
```

## How Config Is Loaded

- **Main** (or test-case setup) typically loads the path given on the command line (e.g. `../json_Files/Test_Config_AMR.json`).
- If the path is a JSON file, **Configuration_Read** (or equivalent) parses it and sets globals from `TestCase`, `Simulation`, `Solver`, and `LimiterCoefficients`.
- AMR keys are **optional**: if absent, defaults are used (`Enable_AMR=false`, `AMR_Period=100`, etc.).
