# Solutions (run output)

Kept production fields and figures. Smoke VTK, duplicate sweep copies, and superseded restarts were removed.

| Path | Role |
|------|------|
| `2D_Euler_Solutions/Cavity_3D/` | Mach-6 cavity LES (block mesh + earlier body-fitted) |
| `2D_Euler_Solutions/SWBLI_3D/` | Laminar plate (`swbli_plate_fine.vtk`) and WALE LES |
| `2D_Euler_Solutions/Ramp_15_Degree_3D/` | Finest 3D ramp 100k |
| `2D_Euler_Solutions/Ramp_15_Degree/` | 2D ramp (fine 481×161 + comparison PNGs) |
| `2D_Euler_Solutions/Half_Cylinder/` | M=6 half-cylinder GridSize_7 |
| `plots/` | Half-cylinder validation figures |
| `les_plots/`, `ramp_plots/`, `ramp_comparison_plots/` | Extra figures |
| `results/` | Selected continue dumps / WENO-order logs (not duplicate VTK) |
