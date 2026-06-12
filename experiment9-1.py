"""
plot_experiment_9.py  —  Experiment 9: Preconditioner comparison

Three compact PNGs, all in the style of plot_experiment_8.py:
  exp9_plot_A.png  — kappa(J) and kappa(M^{-1}J) vs h  (1x2 panels)
  exp9_plot_B.png  — iteration counts (left) and speedup (right), side by side

Note: SSOR is excluded — it requires explicit lower-triangle entries but
the solver stores only the upper triangle (symmetric CSR), so the forward
sweep degenerates to diagonal scaling and SSOR reduces to Jacobi.

Usage:  python plot_experiment_9.py [--csv PATH]
"""

import argparse, csv, math
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--csv", default="experiment_9_preconditioners.csv")
args = parser.parse_args()

rows = []
with open(args.csv, newline="") as fh:
    for r in csv.DictReader(fh):
        rows.append({k: float(v) for k, v in r.items()})

rows.sort(key=lambda r: r["h"], reverse=True)   # coarse -> fine left to right

h         = np.array([r["h"]                  for r in rows])
kappa_J   = np.array([r["kappa_J"]            for r in rows])
kappa_pJ  = np.array([r["kappa_precond_jacobi"] for r in rows])
kappa_pIC = np.array([r["kappa_precond_ic"]     for r in rows])
iters_cg  = np.array([r["iters_cg"]            for r in rows])
iters_J   = np.array([r["iters_jacobi"]        for r in rows])
iters_IC  = np.array([r["iters_ic"]            for r in rows])
spdup_J   = np.array([r["speedup_jacobi"]      for r in rows])
spdup_IC  = np.array([r["speedup_ic"]          for r in rows])

# ---- colour palette (same as experiment 8) --------------------------------
BLUE   = "#2166ac"
RED    = "#d6604d"
ORANGE = "#f4a742"
GREEN  = "#4dac26"
GREY   = "#555555"
BROWN  = "#8b4513"

FS = 12    # axis label font size
TS = 10    # tick label font size
LW = 2.8   # line width
MS = 8     # marker size

# ---------------------------------------------------------------------------
# Plot A — condition numbers vs h
#   Left panel:  kappa(J)          with h^{-2} reference
#   Right panel: kappa(M^{-1}J) for the three preconditioners
# ---------------------------------------------------------------------------
k_ref2 = kappa_J[0] * (h / h[0])**(-2)

fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(9, 3.6))
fig.subplots_adjust(left=0.10, right=0.97, bottom=0.17, top=0.87, wspace=0.38)

# left: kappa(J)
ax0.loglog(h, k_ref2,  color=GREY, ls="--", lw=2.0, label=r"$\propto h^{-2}$")
ax0.loglog(h, kappa_J, color=BLUE, ls="-",  marker="o", ms=MS, lw=LW,
           markerfacecolor="white", markeredgewidth=2.2,
           label=r"$\kappa_2(J)$")
ax0.invert_xaxis()
for i in range(1, len(rows)):
    h_mid = math.sqrt(h[i-1] * h[i])
    k_mid = math.sqrt(kappa_J[i-1] * kappa_J[i])
    slope = math.log(kappa_J[i] / kappa_J[i-1]) / math.log(h[i] / h[i-1])
    ax0.annotate(f"{slope:.2f}", xy=(h_mid, k_mid),
                 xytext=(5, 4), textcoords="offset points",
                 fontsize=TS, color=BLUE)
ax0.set_xlabel(r"Mesh size $h$", fontsize=FS)
ax0.set_ylabel(r"$\kappa_2(J)$", fontsize=FS)
ax0.set_title(r"$\kappa_2(J) \sim h^{-2}$", fontsize=FS)
ax0.tick_params(labelsize=TS)
ax0.legend(fontsize=TS, framealpha=0.92)
ax0.grid(True, which="both", ls="--", lw=0.4, alpha=0.5)

# right: kappa(M^{-1}J) for Jacobi and IC(0)
ax1.loglog(h, kappa_pJ,  color=ORANGE, ls="-",  marker="o", ms=MS, lw=LW,
           markerfacecolor="white", markeredgewidth=2.2, label=r"Jacobi")
ax1.loglog(h, kappa_pIC, color=RED,    ls="-",  marker="^", ms=MS, lw=LW,
           markerfacecolor="white", markeredgewidth=2.2, label=r"IC(0)")
ax1.loglog(h, k_ref2,    color=GREY,   ls="--", lw=2.0,
           label=r"$\propto h^{-2}$")
ax1.invert_xaxis()
for arr, color in [(kappa_pJ, ORANGE), (kappa_pIC, RED)]:
    for i in range(1, len(rows)):
        h_mid = math.sqrt(h[i-1] * h[i])
        k_mid = math.sqrt(arr[i-1] * arr[i])
        slope = math.log(arr[i] / arr[i-1]) / math.log(h[i] / h[i-1])
        ax1.annotate(f"{slope:.2f}", xy=(h_mid, k_mid),
                     xytext=(5, 4), textcoords="offset points",
                     fontsize=TS, color=color)
ax1.set_xlabel(r"Mesh size $h$", fontsize=FS)
ax1.set_ylabel(r"$\kappa_2(M^{-1}J)$", fontsize=FS)
ax1.set_title(r"Preconditioned spectrum", fontsize=FS)
ax1.tick_params(labelsize=TS)
ax1.legend(fontsize=TS, framealpha=0.92)
ax1.grid(True, which="both", ls="--", lw=0.4, alpha=0.5)

fig.savefig("exp9_plot_A.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved exp9_plot_A.png")

# ---------------------------------------------------------------------------
# Plot B — iteration counts (left) and speedup (right), side by side
# ---------------------------------------------------------------------------
cg_ref   = iters_cg[0] * (h / h[0])**(-1)
ic_ref   = iters_IC[0] * (h / h[0])**(-0.5)   # theory: IC(0) gives h^{-1/2}
ic_slope = np.polyfit(np.log(h), np.log(iters_IC), 1)[0]  # observed slope for annotations

fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.0))
fig.subplots_adjust(left=0.09, right=0.97, bottom=0.14, top=0.75, wspace=0.36)

# --- left: iteration counts ---
ax0.loglog(h, cg_ref,   color=GREY,   ls="--", lw=2.0, label=r"$\propto h^{-1}$")
ax0.loglog(h, ic_ref,   color=GREY,   ls=":",  lw=2.0, label=r"$\propto h^{-1/2}$ (theory)")
ax0.loglog(h, iters_cg, color=BLUE,   ls="-",  marker="o", ms=MS,   lw=LW,
           markerfacecolor="white", markeredgewidth=2.2, label="CG (no precond)")
ax0.loglog(h, iters_J,  color=ORANGE, ls="-",  marker="s", ms=MS,   lw=LW,
           markerfacecolor="white", markeredgewidth=2.2, label="PCG-Jacobi")
ax0.loglog(h, iters_IC, color=RED,    ls="-",  marker="D", ms=MS-1, lw=LW,
           markerfacecolor="white", markeredgewidth=2.2, label="PCG-IC(0)")

for arr, color in [(iters_cg, BLUE), (iters_IC, RED)]:
    for i in range(1, len(rows)):
        if arr[i] > 0 and arr[i-1] > 0:
            h_mid  = math.sqrt(h[i-1] * h[i])
            it_mid = math.sqrt(arr[i-1] * arr[i])
            slope  = math.log(arr[i] / arr[i-1]) / math.log(h[i] / h[i-1])
            ax0.annotate(f"{slope:.2f}", xy=(h_mid, it_mid),
                         xytext=(5, 4), textcoords="offset points",
                         fontsize=TS - 1, color=color)

ax0.invert_xaxis()
ax0.set_xlabel(r"Mesh size $h$", fontsize=FS)
ax0.set_ylabel("Iterations", fontsize=FS)
ax0.set_title("Iteration count vs mesh size", fontsize=FS)
ax0.tick_params(labelsize=TS)
ax0.legend(fontsize=TS, framealpha=0.92, loc="upper center",
           bbox_to_anchor=(0.5, 1.42), ncol=2)
ax0.grid(True, which="both", ls="--", lw=0.4, alpha=0.5)

# --- right: speedup ---
ax1.semilogx(h, spdup_J,  color=ORANGE, ls="-",  marker="o", ms=MS, lw=LW,
             markerfacecolor="white", markeredgewidth=2.2, label="PCG-Jacobi")
ax1.semilogx(h, spdup_IC, color=RED,    ls="-",  marker="^", ms=MS, lw=LW,
             markerfacecolor="white", markeredgewidth=2.2, label="PCG-IC(0)")
ax1.axhline(1.0, color=GREY, ls="--", lw=1.5, label="no speedup")

ax1.invert_xaxis()
ax1.set_xlabel(r"Mesh size $h$", fontsize=FS)
ax1.set_ylabel(r"Speedup $= n_{\mathrm{CG}} / n_{\mathrm{PCG}}$", fontsize=FS)
ax1.set_title("Speedup over unpreconditioned CG", fontsize=FS)
ax1.tick_params(labelsize=TS)
ax1.legend(fontsize=TS, framealpha=0.92, loc="upper center",
           bbox_to_anchor=(0.5, 1.42), ncol=2)
ax1.grid(True, ls="--", lw=0.4, alpha=0.5)

fig.savefig("exp9_plot_B.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved exp9_plot_B.png")