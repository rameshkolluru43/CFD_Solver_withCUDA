# Half-Cylinder M∞ = 6 Validation Guide

Mach-6 flow over a half cylinder on a structured 481×161 (480×160 cell) grid. Use this page for the canonical configs, expected numbers, and plotting after the July 2026 viscous corner / face-normal fix.

See also: [overview.md](../overview.md), [MESH_AND_GRID.md](MESH_AND_GRID.md), [RELEASE_NOTES.md](RELEASE_NOTES.md).

## Grid and conventions

| Item | Value |
|------|--------|
| Grid file | `Grid_Files/Half_Cylinder/Half_Cylinder_481_161.txt` |
| Cells | 480 × 160 |
| Face order | **LBRT** — neighbours `[left, bottom, right, top]` |
| Vertices | File order `[o, a, b, c]` (do **not** atan2-reorder) |

Wall BC faces must keep `Face_Normals` aligned with Neighbours / boundary lists. Misalignment at the side corners previously produced viscous blow-up.

## Case and run JSON

| File | Role |
|------|------|
| `input/json_Files/Half_Cylinder_481_161_M6.json` | Inviscid freestream case params |
| `input/json_Files/Half_Cylinder_481_161_M6_freestreamIC.json` | Freestream IC variant |
| `input/json_Files/Half_Cylinder_481_161_M6_viscous.json` | Viscous: isothermal no-slip, Re = 1e5, Tw = 2 |
| `input/json_Files/Run_HalfCylinder_M6_corner_fix_smoke500.json` | P3 → viscous 500-iter smoke (`Initialize_Type=1`) |
| `input/json_Files/Run_HalfCylinder_M6_WENO_RICCA_P3_viscous_20k.json` | Longer viscous continue |
| `input/json_Files/Run_HalfCylinder_M6_WENO_RICCA_P3_viscous_smoke.json` | Short viscous smoke |
| `input/json_Files/Run_HalfCylinder_M6_WENO_RICCA_inviscid_smoke_host.json` | Host inviscid WENO+RICCA smoke |

Typical solver settings for the validated path:

- `Is_WENO: true`, `Dissipation_Type: 4` (RICCA)
- Viscous: `Is_Viscous: true` (sets viscous wall / NS path)
- Restart smoke: `Initialize_Type: 1`, `CFL: 0.02`

## Build and run

```bash
cmake -S . -B build-cuda
cmake --build build-cuda --target CFD_solver_gpu -j$(nproc)

# Viscous corner-fix smoke (continue from P3 solution)
./build-cuda/CFD_solver_gpu input/json_Files/Run_HalfCylinder_M6_corner_fix_smoke500.json
```

Host CPU build (`CFD_solver`) is also valid for the WENO+RICCA path; CUDA is optional for this validation.

## Expected results (reference)

| Case | Location | Pmax | Mmax |
|------|----------|------|------|
| Inviscid WENO+RICCA | `results/ordered_sweep/P3_RICCA_WENO/` | ≈ 48.3 | ≈ 6 |
| Viscous smoke 500 iters after corner fix | `results/halfcyl_corner_fix_smoke500/` | ≈ 48.9 | ≈ 6.05 |
| Pre-fix viscous continue (broken) | — | ~137 | ~64 |

## Plots and diagnostics

```bash
python3 scripts/plot_halfcyl_inviscid_viscous.py
```

Outputs under `plots/half_cylinder_inviscid_viscous/`:

| Artifact | Path |
|----------|------|
| Mach / pressure / residuals | `*_mach.png`, `*_pressure.png`, `*_residuals.png`, `*_summary.png` |
| Side-by-side compare | `compare_mach_pressure.png`, `compare_residuals.png` |
| Corner face/normal dump | `viscous_corner_diagnostics.txt` |
| Status note | `README_status.txt` |

When viscous, the solver also writes `viscous_corner_diagnostics.txt` under the solution directory via `Dump_Viscous_Corner_Diagnostics`.

## Related robustness fixes (same campaign)

1. Structured TXT: keep `[o,a,b,c]`; LBRT face normals outward-corrected.
2. Corner secondary neighbours for co-volume gradients.
3. Supersonic inlet: prescribe freestream (do not mutate `inletCond`).
4. Host WENO wall faces: 1O fallback on degenerate stencils.
5. `Update()`: reject non-positive updates instead of flooring pressure.
6. Viscous driver: remove double first-order flux evaluation before the explicit NS step.

## WENO order-of-accuracy helper

```bash
python3 scripts/test_weno_order_accuracy.py
```

Writes under `plots/weno_order_accuracy/`. Complementary to the half-cylinder campaign; see [WENO2D_Validation_Report.md](WENO2D_Validation_Report.md).
