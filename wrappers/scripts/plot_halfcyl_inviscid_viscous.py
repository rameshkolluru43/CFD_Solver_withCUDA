#!/usr/bin/env python3
"""Plot inviscid vs viscous half-cylinder WENO + RICCA (post corner-normal fix)."""
from __future__ import annotations

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

ROOT = Path("/home/ramesh.kolluru/CFD_Solver_withCUDA")
PLOT = ROOT / "plots" / "half_cylinder_inviscid_viscous"
NX, NY = 481, 161
NCX, NCY = NX - 1, NY - 1
NPTS, NCELLS = NX * NY, NCX * NCY
GAMMA = 1.4

CASES = [
    (
        "Inviscid WENO + RICCA (P3)",
        ROOT / "results/ordered_sweep/P3_RICCA_WENO",
        ROOT / "results/ordered_sweep/P3_RICCA_WENO/tail.log",
        "inviscid",
        "P3 reference",
    ),
    (
        "Viscous WENO + RICCA + NS (corner fix)",
        ROOT / "results/halfcyl_corner_fix_smoke500",
        ROOT / "results/halfcyl_corner_fix_smoke500.log",
        "viscous",
        "500 iters from P3",
    ),
]


def section_start(lines, tag):
    for i, ln in enumerate(lines):
        toks = ln.replace("\t", " ").split()
        if toks and toks[0] == tag:
            return i, toks
    raise RuntimeError(f"missing section {tag}")


def read_floats_after(lines, start_idx, nvals):
    vals = []
    i = start_idx
    while len(vals) < nvals:
        for p in lines[i].replace("\t", " ").split():
            try:
                vals.append(float(p))
            except ValueError:
                continue
            if len(vals) >= nvals:
                break
        i += 1
    return np.asarray(vals[:nvals], dtype=float), i


def read_vtk(path: Path):
    lines = path.read_text().splitlines()
    i_pts, toks = section_start(lines, "POINTS")
    npts = int(toks[1])
    pts, _ = read_floats_after(lines, i_pts + 1, npts * 3)
    pts = pts.reshape(npts, 3)

    x = pts[:NPTS, 0].reshape(NY, NX)
    y = pts[:NPTS, 1].reshape(NY, NX)

    fields = {}
    i = 0
    while i < len(lines):
        toks = lines[i].replace("\t", " ").split()
        if toks and toks[0] == "SCALARS" and len(toks) >= 2:
            name = toks[1]
            j = i + 1
            while j < len(lines) and not lines[j].replace("\t", " ").split():
                j += 1
            if j < len(lines) and lines[j].replace("\t", " ").split()[0] == "LOOKUP_TABLE":
                j += 1
            vals, _ = read_floats_after(lines, j, NCELLS)
            if name not in ("AMR_Level", "AMR_Parent", "AMR_IsLeaf"):
                fields[name] = vals.reshape(NCY, NCX)
            i = j
        i += 1
    return x, y, fields


def parse_residuals(log_path: Path):
    iters, rho, rhou, rhov, rhoe = [], [], [], [], []
    if not log_path.exists():
        return None
    for ln in log_path.read_text().splitlines():
        parts = ln.split()
        if len(parts) < 6:
            continue
        try:
            it = int(parts[0])
            if it <= 0:
                continue
            e0, e1, e2, e3 = map(float, parts[2:6])
        except ValueError:
            continue
        iters.append(it)
        rho.append(e0)
        rhou.append(e1)
        rhov.append(e2)
        rhoe.append(e3)
    if not iters:
        return None
    return {
        "iter": np.asarray(iters),
        "rho": np.asarray(rho),
        "rhou": np.asarray(rhou),
        "rhov": np.asarray(rhov),
        "rhoe": np.asarray(rhoe),
    }


def plot_field(ax, x, y, z, title, cmap="turbo", vmin=None, vmax=None):
    if vmin is None:
        vmin = np.nanpercentile(z, 1)
    if vmax is None:
        vmax = np.nanpercentile(z, 99)
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin >= vmax:
        vmin, vmax = float(np.nanmin(z)), float(np.nanmax(z))
        if vmin >= vmax:
            vmax = vmin + 1e-12
    pcm = ax.pcolormesh(
        x, y, z, shading="flat", cmap=cmap, norm=Normalize(vmin=vmin, vmax=vmax)
    )
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    return pcm


def field_stats(fields):
    p = fields["Pressure"]
    m = fields["Mach_Number"]
    return float(np.nanmax(p)), float(np.nanmax(m))


def main():
    PLOT.mkdir(parents=True, exist_ok=True)
    for p in PLOT.glob("*.png"):
        p.unlink()

    loaded = []
    for label, case_dir, log_path, tag, note in CASES:
        vtk = case_dir / "Half_Cylinder_Grid_Size_7_FluxScheme_WENO_Dissipation_RICCA.vtk"
        if not vtk.exists():
            raise FileNotFoundError(vtk)
        xc, yc, fields = read_vtk(vtk)
        if "Mach_Number" not in fields and all(
            k in fields for k in ("Density", "Pressure", "U_Velocity", "V_Velocity")
        ):
            a = np.sqrt(
                np.maximum(GAMMA * fields["Pressure"] / np.maximum(fields["Density"], 1e-14), 0.0)
            )
            fields["Mach_Number"] = np.hypot(fields["U_Velocity"], fields["V_Velocity"]) / np.maximum(
                a, 1e-14
            )
        pmax, mmax = field_stats(fields)
        res = parse_residuals(log_path)
        loaded.append((label, tag, note, xc, yc, fields, res, case_dir, pmax, mmax))

        fig, axs = plt.subplots(2, 2, figsize=(10.5, 8.5), constrained_layout=True)
        keys = [
            ("Mach_Number", "Mach", "turbo", 0.0, 6.0),
            ("Pressure", "Pressure", "viridis", None, None),
            ("Density", "Density", "plasma", None, None),
            ("Temperature", "Temperature", "inferno", None, None),
        ]
        for ax, (key, title, cmap, vmin, vmax) in zip(axs.ravel(), keys):
            if key not in fields:
                ax.set_visible(False)
                continue
            pcm = plot_field(ax, xc, yc, fields[key], title, cmap=cmap, vmin=vmin, vmax=vmax)
            fig.colorbar(pcm, ax=ax, fraction=0.046, pad=0.04)
        fig.suptitle(f"{label} — {note}\nPmax={pmax:.2f}, Mmax={mmax:.2f}", fontsize=12)
        fig.savefig(PLOT / f"{tag}_summary.png", dpi=160)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(6.2, 5.2), constrained_layout=True)
        pcm = plot_field(
            ax, xc, yc, fields["Mach_Number"], f"{label}\nMach (Mmax={mmax:.2f})", cmap="turbo", vmin=0, vmax=6
        )
        fig.colorbar(pcm, ax=ax, fraction=0.046, pad=0.04, label="M")
        fig.savefig(PLOT / f"{tag}_mach.png", dpi=160)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(6.2, 5.2), constrained_layout=True)
        pcm = plot_field(
            ax, xc, yc, fields["Pressure"], f"{label}\nPressure (Pmax={pmax:.2f})", cmap="viridis"
        )
        fig.colorbar(pcm, ax=ax, fraction=0.046, pad=0.04, label="P")
        fig.savefig(PLOT / f"{tag}_pressure.png", dpi=160)
        plt.close(fig)

        if res is not None:
            fig, ax = plt.subplots(figsize=(7.2, 4.2), constrained_layout=True)
            ax.semilogy(res["iter"], res["rho"], label=r"$\rho$")
            ax.semilogy(res["iter"], res["rhov"], label=r"$\rho v$")
            ax.semilogy(res["iter"], res["rhoe"], label=r"$\rho E$")
            ax.semilogy(res["iter"], res["rhou"], label=r"$\rho u$", alpha=0.45)
            ax.set_xlabel("Iteration")
            ax.set_ylabel("Residual")
            ax.set_title(f"{label} residuals")
            ax.grid(True, which="both", alpha=0.3)
            ax.legend(fontsize=8)
            fig.savefig(PLOT / f"{tag}_residuals.png", dpi=160)
            plt.close(fig)

    fig, axs = plt.subplots(2, 2, figsize=(11, 9), constrained_layout=True)
    for col, (label, tag, note, xc, yc, fields, res, case_dir, pmax, mmax) in enumerate(loaded):
        pcm0 = plot_field(
            axs[0, col],
            xc,
            yc,
            fields["Mach_Number"],
            f"{tag}\nMach (Mmax={mmax:.2f})",
            cmap="turbo",
            vmin=0,
            vmax=6,
        )
        fig.colorbar(pcm0, ax=axs[0, col], fraction=0.046, pad=0.04)
        pcm1 = plot_field(
            axs[1, col],
            xc,
            yc,
            fields["Pressure"],
            f"{tag}\nPressure (Pmax={pmax:.2f})",
            cmap="viridis",
        )
        fig.colorbar(pcm1, ax=axs[1, col], fraction=0.046, pad=0.04)
    fig.suptitle(
        "Half-cylinder M∞=6 — inviscid P3 vs viscous (corner face/normal fix)", fontsize=12
    )
    fig.savefig(PLOT / "compare_mach_pressure.png", dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.5, 4.4), constrained_layout=True)
    for label, tag, note, xc, yc, fields, res, case_dir, pmax, mmax in loaded:
        if res is None:
            continue
        ax.semilogy(res["iter"], res["rho"], label=f"{tag} ρ")
        ax.semilogy(res["iter"], res["rhov"], "--", label=f"{tag} ρv")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Residual")
    ax.set_title("Residual history")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    fig.savefig(PLOT / "compare_residuals.png", dpi=160)
    plt.close(fig)

    (PLOT / "README_status.txt").write_text(
        "Replotted 2026-07-20 after corner face/normal fix.\n"
        "Inviscid: results/ordered_sweep/P3_RICCA_WENO\n"
        "Viscous:  results/halfcyl_corner_fix_smoke500 (500 iters from P3)\n"
        "Removed prior bad viscous_mach / P3_viscous_20k_vs_P3 plots.\n"
    )

    print("Wrote plots to", PLOT)
    for p in sorted(PLOT.glob("*.png")):
        print(" ", p.name, f"({p.stat().st_size} bytes)")
    for label, tag, note, xc, yc, fields, res, case_dir, pmax, mmax in loaded:
        print(f"{tag}: {note}  Pmax={pmax:.3f}  Mmax={mmax:.3f}")


if __name__ == "__main__":
    main()
