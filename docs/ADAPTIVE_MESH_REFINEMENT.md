# Adaptive Mesh Refinement (AMR)

The solver includes **gradient-based dynamic meshing**: cells are tagged for refinement using solution gradients. The current implementation performs **tagging only** (no mesh topology change).

## Overview

- **Indicator**: Per-cell scalar combining density and pressure gradients, made scale-invariant with cell size.
- **Tagging**: Cells with indicator above a threshold are marked `Is_Splittable`; optionally only a fraction of cells (by indicator rank) are tagged.
- **Solver hook**: Every **AMR_Period** iterations, the indicator is computed and cells are tagged; a short log line reports how many cells were tagged.

## Gradient-Based Indicator

For each cell:

1. **Green–Gauss gradient** at the cell centre (any polygon):
   - \( \nabla\phi = \frac{1}{A} \sum_f \bar\phi_f \, \mathbf{n}_f \, \Delta\ell_f \)
   - \( \bar\phi_f = \frac{1}{2}(\phi_{\text{cell}} + \phi_{\text{neighbour}}) \), or \( \phi_{\text{cell}} \) on boundary faces.

2. **Refinement indicator** (stored in `Gradient_Refinement_Indicator[i]`):
   - \( h = \sqrt{A} \)
   - \( \eta = |\nabla\rho|\, h + \frac{|\nabla P|}{\max(P,\epsilon)}\, h \)

So the indicator is large where density or pressure change quickly relative to cell size.

## Tagging Rules

- **Threshold**: `AMR_Gradient_Threshold`. Cells with \( \eta > \text{threshold} \) are candidates.
- **Cap (optional)**: If `AMR_Max_Fraction` is in \( (0,1) \), only the **top AMR_Max_Fraction** of cells (by \( \eta \)) are tagged. This limits the number of refined cells per pass.

After tagging, `Cells[i].Is_Splittable` is `true` for tagged cells.

## Configuration (JSON)

Under the `Solver` section:

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `Enable_AMR` | bool | false | Turn on gradient-based AMR tagging. |
| `AMR_Period` | int | 100 | Run AMR (compute indicator + tag) every this many iterations. |
| `AMR_Gradient_Threshold` | double | 0.1 | Refine cells with indicator above this value. |
| `AMR_Max_Fraction` | double | 0.3 | If in (0,1), refine at most this fraction of cells (by indicator). 0 or ≥1 = no cap. |

Example (see `json_Files/Test_Config_AMR.json`):

```json
"Solver": {
  "Enable_AMR": true,
  "AMR_Period": 100,
  "AMR_Gradient_Threshold": 0.05,
  "AMR_Max_Fraction": 0.25
}
```

## Implementation Details

- **Compute_Gradient_Refinement_Indicator()**: Fills `Gradient_Refinement_Indicator` for all physical cells.
- **TagRefinableCells(cells, threshold)**: Uses the indicator and threshold (and optional max fraction) to set `Is_Splittable`.
- **Apply_Adaptive_Refinement()**: Called from the inviscid solver every `AMR_Period`; it calls `TagRefinableCells` and logs the number of tagged cells. It returns `false` (no mesh change).

## Future: Actual Refinement

Planned extension:

- For cells with `Is_Splittable` (e.g. quads), **split** into 4 children (edge midpoints + centre).
- Build new global point list and new cell list; update connectivity and boundary lists.
- **Solution transfer**: copy or interpolate conservative state from parent to children; resize state arrays and re-run initialization/ghost construction.

Until then, AMR is **indicator + tagging only**; the tagged cells can be used for visualization or for an external refinement tool.
