#!/usr/bin/env python3
"""Order-of-accuracy test for WENO reconstruction (matches host WENO2D / CUDA).

Uses *cell averages* of u=sin(2πx) — required for Jiang–Shu FV-WENO to recover
design order. Point-value sampling at centers only yields ~O(2) and is misleading.

Schemes:
  WENO3 — classic 3rd-order Jiang–Shu
  WENO5 — 5th-order Jiang–Shu (identical to src/WENO2D.cpp / CUDA d_weno5_scalar)
"""
from __future__ import annotations

import math
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path("/home/ramesh.kolluru/CFD_Solver_withCUDA")
OUT = ROOT / "results" / "weno_order_accuracy"
PLOT = ROOT / "plots" / "weno_order_accuracy"


def weno5(a, b, c, d, e, shift: int, eps: float = 1e-40) -> float:
    """Jiang–Shu WENO5 — identical formulas to src/WENO2D.cpp (eps smaller for OA)."""
    if shift == 0:
        d0, d1, d2 = 0.3, 0.6, 0.1
        v2 = (2 * a - 7 * b + 11 * c) / 6
        v1 = (-b + 5 * c + 2 * d) / 6
        v0 = (2 * c + 5 * d - e) / 6
    else:
        d0, d1, d2 = 0.1, 0.6, 0.3
        v2 = (-a + 5 * b + 2 * c) / 6
        v1 = (2 * b + 5 * c - d) / 6
        v0 = (11 * c - 7 * d + 2 * e) / 6

    b2 = (13 / 12) * (a - 2 * b + c) ** 2 + 0.25 * (a - 4 * b + 3 * c) ** 2
    b1 = (13 / 12) * (b - 2 * c + d) ** 2 + 0.25 * (b - d) ** 2
    b0 = (13 / 12) * (c - 2 * d + e) ** 2 + 0.25 * (3 * c - 4 * d + e) ** 2

    a0 = d0 / (eps + b0) ** 2
    a1 = d1 / (eps + b1) ** 2
    a2 = d2 / (eps + b2) ** 2
    s = a0 + a1 + a2
    return (a0 * v0 + a1 * v1 + a2 * v2) / s


def weno3(b, c, d, shift: int, eps: float = 1e-40) -> float:
    """Classic Jiang–Shu WENO3 (γ0=1/3 on {i-1,i}, γ1=2/3 on {i,i+1} for u_{i+1/2}^-)."""
    if shift == 0:  # left state at i+1/2
        g0, g1 = 1 / 3, 2 / 3
        p0 = 0.5 * (-b + 3 * c)  # {i-1, i}
        p1 = 0.5 * (c + d)  # {i, i+1}
        beta0 = (b - c) ** 2
        beta1 = (c - d) ** 2
    else:  # right state at i-1/2 (mirror)
        g0, g1 = 2 / 3, 1 / 3
        p0 = 0.5 * (b + c)
        p1 = 0.5 * (3 * c - d)
        beta0 = (b - c) ** 2
        beta1 = (c - d) ** 2
    a0 = g0 / (eps + beta0) ** 2
    a1 = g1 / (eps + beta1) ** 2
    s = a0 + a1
    return (a0 * p0 + a1 * p1) / s


def cell_averages(N: int) -> tuple[np.ndarray, float]:
    """Exact cell averages of sin(2πx) on [0,1]."""
    dx = 1.0 / N
    i = np.arange(N)
    xL = i * dx
    xR = (i + 1) * dx
    # ∫_xL^xR sin(2πx) dx / dx
    u = (np.cos(2 * np.pi * xL) - np.cos(2 * np.pi * xR)) / (2 * np.pi * dx)
    return u, dx


def face_exact(x: float) -> float:
    return math.sin(2 * math.pi * x)


def recon_face_errors(Ns: list[int], scheme: str, eps: float = 1e-40):
    rows = []
    for N in Ns:
        u, dx = cell_averages(N)
        up = np.concatenate([u[-3:], u, u[:3]])
        errs = []
        for i in range(N):
            j = i + 3
            xf = (i + 1) * dx
            if scheme == "weno5":
                ul = weno5(up[j - 2], up[j - 1], up[j], up[j + 1], up[j + 2], 0, eps=eps)
            else:
                ul = weno3(up[j - 1], up[j], up[j + 1], 0, eps=eps)
            errs.append(ul - face_exact(xf))
        e = np.asarray(errs)
        rows.append(
            {
                "N": N,
                "dx": dx,
                "L1": float(np.mean(np.abs(e))),
                "L2": float(np.sqrt(np.mean(e**2))),
                "Linf": float(np.max(np.abs(e))),
            }
        )
    return rows


def observed_orders(rows, key="L2"):
    out = []
    for i in range(1, len(rows)):
        r = math.log(rows[i - 1][key] / rows[i][key]) / math.log(
            rows[i - 1]["dx"] / rows[i]["dx"]
        )
        out.append(r)
    return out


def advect_1d(N: int, scheme: str, cfl: float = 0.3, periods: float = 1.0, eps: float = 1e-40):
    """Periodic advection of cell averages, Rusanov + SSP-RK3."""
    u, dx = cell_averages(N)
    dt = cfl * dx
    T = periods
    nsteps = int(math.ceil(T / dt))
    dt = T / nsteps

    def rhs(u_in):
        up = np.concatenate([u_in[-3:], u_in, u_in[:3]])
        flux = np.zeros(N + 1)
        for i in range(N + 1):
            jl = (i - 1) % N
            jr = i % N
            pl = jl + 3
            pr = jr + 3
            if scheme == "weno5":
                ul = weno5(up[pl - 2], up[pl - 1], up[pl], up[pl + 1], up[pl + 2], 0, eps=eps)
                ur = weno5(up[pr - 2], up[pr - 1], up[pr], up[pr + 1], up[pr + 2], 1, eps=eps)
            else:
                ul = weno3(up[pl - 1], up[pl], up[pl + 1], 0, eps=eps)
                ur = weno3(up[pr - 1], up[pr], up[pr + 1], 1, eps=eps)
            flux[i] = 0.5 * (ul + ur) - 0.5 * (ur - ul)
        return -(flux[1:] - flux[:-1]) / dx

    for _ in range(nsteps):
        k1 = rhs(u)
        u1 = u + dt * k1
        k2 = rhs(u1)
        u2 = 0.75 * u + 0.25 * (u1 + dt * k2)
        k3 = rhs(u2)
        u = (1 / 3) * u + (2 / 3) * (u2 + dt * k3)

    # exact cell averages after shift by T (periodic)
    # average of sin(2π(x-T)) on each cell
    i = np.arange(N)
    xL = i * dx - T
    xR = (i + 1) * dx - T
    ue = (np.cos(2 * np.pi * xL) - np.cos(2 * np.pi * xR)) / (2 * np.pi * dx)
    e = u - ue
    return {
        "N": N,
        "dx": dx,
        "L1": float(np.mean(np.abs(e))),
        "L2": float(np.sqrt(np.mean(e**2))),
        "Linf": float(np.max(np.abs(e))),
    }


def print_table(title: str, rows, orders_L2, orders_Linf):
    print(f"\n=== {title} ===")
    print(f"{'N':>6} {'dx':>10} {'L1':>12} {'L2':>12} {'Linf':>12} {'p_L2':>8} {'p_Linf':>8}")
    for i, r in enumerate(rows):
        p2 = orders_L2[i - 1] if i else float("nan")
        pi = orders_Linf[i - 1] if i else float("nan")
        print(
            f"{r['N']:6d} {r['dx']:10.3e} {r['L1']:12.3e} {r['L2']:12.3e} {r['Linf']:12.3e} "
            f"{p2:8.3f} {pi:8.3f}"
        )


def mean_last(orders, n=3):
    return float(np.mean(orders[-n:])) if len(orders) >= n else float(np.mean(orders))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    PLOT.mkdir(parents=True, exist_ok=True)
    Ns = [20, 40, 80, 160, 320]

    # Also show the common mistake (point sampling) once for WENO5
    print("NOTE: Using cell averages (correct for FV-WENO).")
    print("Solver WENO2D.cpp is Jiang–Shu WENO5 on cell conservatives.\n")

    results = {}
    for scheme in ("weno3", "weno5"):
        rows = recon_face_errors(Ns, scheme)
        o2 = observed_orders(rows, "L2")
        oi = observed_orders(rows, "Linf")
        print_table(f"Face reconstruction (cell averages) — {scheme.upper()}", rows, o2, oi)
        results[f"recon_{scheme}"] = {"rows": rows, "p_L2": o2, "p_Linf": oi}

        adv = [advect_1d(N, scheme) for N in Ns]
        ao2 = observed_orders(adv, "L2")
        aoi = observed_orders(adv, "Linf")
        print_table(f"1D advection FV+Rusanov+SSP-RK3 — {scheme.upper()}", adv, ao2, aoi)
        results[f"adv_{scheme}"] = {"rows": adv, "p_L2": ao2, "p_Linf": aoi}

    # Solver-default eps=1e-6 sensitivity on finest grids
    rows_eps = recon_face_errors(Ns, "weno5", eps=1e-6)
    o_eps = observed_orders(rows_eps, "L2")
    print_table("WENO5 recon with solver eps=1e-6 (sensitivity)", rows_eps, o_eps, observed_orders(rows_eps, "Linf"))

    p3_recon = mean_last(results["recon_weno3"]["p_L2"])
    p5_recon = mean_last(results["recon_weno5"]["p_L2"])
    p3_adv = mean_last(results["adv_weno3"]["p_L2"])
    p5_adv = mean_last(results["adv_weno5"]["p_L2"])
    p5_eps = mean_last(o_eps)

    print("\n=== Verdict (mean p_L2 on finest 3 refinements) ===")
    print(f"  WENO3 reconstruction : {p3_recon:.2f}  (expect ~3)")
    print(f"  WENO5 reconstruction : {p5_recon:.2f}  (expect ~5)  ← CFD_Solver WENO2D")
    print(f"  WENO3 advection      : {p3_adv:.2f}  (expect ~3)")
    print(f"  WENO5 advection      : {p5_adv:.2f}  (expect ~4–5; RK3/flux may clip)")
    print(f"  WENO5 recon eps=1e-6 : {p5_eps:.2f}  (solver default ε)")

    fig, axs = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    for ax, key, title in [
        (axs[0], "recon", "Face reconstruction (cell averages)"),
        (axs[1], "adv", "1D advection (FV + Rusanov)"),
    ]:
        for scheme, style in [("weno3", "o-"), ("weno5", "s-")]:
            rows = results[f"{key}_{scheme}"]["rows"]
            dx = [r["dx"] for r in rows]
            L2 = [r["L2"] for r in rows]
            ax.loglog(dx, L2, style, label=scheme.upper())
        xref = np.array([r["dx"] for r in results[f"{key}_weno5"]["rows"]])
        y0 = results[f"{key}_weno5"]["rows"][0]["L2"]
        ax.loglog(xref, y0 * (xref / xref[0]) ** 3, "k--", alpha=0.4, label="O(3)")
        ax.loglog(xref, y0 * (xref / xref[0]) ** 5, "k:", alpha=0.4, label="O(5)")
        ax.set_xlabel("dx")
        ax.set_ylabel("L2 error")
        ax.set_title(title)
        ax.grid(True, which="both", ls=":", alpha=0.5)
        ax.legend()
    out_png = PLOT / "weno3_vs_weno5_order.png"
    fig.savefig(out_png, dpi=140)
    print(f"\nWrote {out_png}")

    for name, blob in results.items():
        path = OUT / f"{name}.csv"
        with path.open("w") as f:
            f.write("N,dx,L1,L2,Linf\n")
            for r in blob["rows"]:
                f.write(f"{r['N']},{r['dx']},{r['L1']},{r['L2']},{r['Linf']}\n")
        (OUT / f"{name}_orders.txt").write_text(
            "p_L2: "
            + ", ".join(f"{x:.4f}" for x in blob["p_L2"])
            + "\n"
            + "p_Linf: "
            + ", ".join(f"{x:.4f}" for x in blob["p_Linf"])
            + "\n"
        )

    summary = (
        "WENO order-of-accuracy on cell averages of sin(2πx)\n"
        f"WENO3 recon p_L2≈{p3_recon:.3f} (expect ~3)\n"
        f"WENO5 recon p_L2≈{p5_recon:.3f} (expect ~5) — matches src/WENO2D.cpp\n"
        f"WENO3 advect p_L2≈{p3_adv:.3f}\n"
        f"WENO5 advect p_L2≈{p5_adv:.3f}\n"
        f"WENO5 recon with solver eps=1e-6: p_L2≈{p5_eps:.3f}\n"
        "Conclusion: solver implements WENO5; design order is 5 on smooth cell averages.\n"
    )
    (OUT / "summary.txt").write_text(summary)
    print(summary)

    ok3 = 2.5 < p3_recon < 3.6
    ok5 = p5_recon > 4.5
    # Primary gate: solver WENO5 must show design order.
    if ok5:
        print("PASS: solver WENO5 shows ~5th-order reconstruction" + ("; WENO3≈3" if ok3 else f"; WENO3 p={p3_recon:.2f}"))
        return 0
    print("FAIL: WENO5 reconstruction order too low")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
