"""
plot_experiment_10.py  —  Experiment 10: Picard vs Newton convergence rate

Two PNGs:
  exp10_plot_A.png  — residual norm vs iteration for all init conditions,
                      Newton and Picard overlaid (one panel per init)
  exp10_plot_B.png  — lin_ratio (Picard) and quad_ratio (Newton) vs iteration
                      for the warm starts, showing the rate signature

Usage:  python plot_experiment_10.py [--csv PATH]
"""
import argparse, csv, math
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--csv", default="experiment_10_picard_vs_newton.csv")
args = parser.parse_args()

rows = []
with open(args.csv, newline="") as fh:
    for r in csv.DictReader(fh):
        rows.append({
            "solver": r["solver"],
            "init":   r["init"],
            **{k: float(v) for k, v in r.items() if k not in ("solver","init")}
        })

# ---- colour / style palette -------------------------------------------
BLUE   = "#2166ac"
RED    = "#d6604d"
ORANGE = "#f4a742"
GREEN  = "#4dac26"
GREY   = "#555555"

SOLVER_COLOR = {"newton": BLUE,  "picard": RED}
SOLVER_LS    = {"newton": "-",   "picard": "--"}
SOLVER_LABEL = {"newton": "Newton", "picard": "Picard"}

INIT_ORDER = ["cold", "warm_1e-01", "warm_1e-02", "warm_1e-04"]
INIT_LABEL = {
    "cold":       "Cold start",
    "warm_1e-01": r"$\varepsilon = 10^{-1}$",
    "warm_1e-02": r"$\varepsilon = 10^{-2}$",
    "warm_1e-04": r"$\varepsilon = 10^{-4}$",
}

FS = 12
TS = 10
LW = 2.4
MS = 7

def get(data, solver, init):
    return [r for r in data if r["solver"] == solver and r["init"] == init]

# ---------------------------------------------------------------------------
# Plot A — residual norm vs iteration, 4 panels (one per init condition)
# ---------------------------------------------------------------------------
inits_present = []
for init in INIT_ORDER:
    if any(r["init"] == init for r in rows):
        inits_present.append(init)

ncols = 2
nrows = math.ceil(len(inits_present) / ncols)
fig, axes = plt.subplots(nrows, ncols,
                         figsize=(10, 3.6 * nrows),
                         squeeze=False)
fig.subplots_adjust(left=0.09, right=0.97, bottom=0.12,
                    top=0.88, wspace=0.32, hspace=0.48)

for idx, init in enumerate(inits_present):
    ax = axes[idx // ncols][idx % ncols]
    for solver in ["newton", "picard"]:
        d = get(rows, solver, init)
        if not d:
            continue
        iters = [r["iter"] for r in d]
        res   = [r["residual_norm"] for r in d]
        ax.semilogy(iters, res,
                    color=SOLVER_COLOR[solver],
                    ls=SOLVER_LS[solver], lw=LW,
                    marker="o" if solver == "newton" else "s",
                    ms=MS, markerfacecolor="white", markeredgewidth=2.0,
                    label=SOLVER_LABEL[solver])
    ax.set_title(INIT_LABEL.get(init, init), fontsize=FS)
    ax.set_xlabel("Iteration $k$", fontsize=FS)
    ax.set_ylabel(r"$\|F(u^k)\|_2$", fontsize=FS)
    ax.tick_params(labelsize=TS)
    ax.legend(fontsize=TS, framealpha=0.92)
    ax.grid(True, which="both", ls="--", lw=0.4, alpha=0.5)

# hide unused panels
for idx in range(len(inits_present), nrows * ncols):
    axes[idx // ncols][idx % ncols].set_visible(False)

fig.savefig("exp10_plot_A.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved exp10_plot_A.png")

# ---------------------------------------------------------------------------
# Plot B — convergence rate diagnostics for warm starts
#   Left:  lin_ratio  for Picard  (should converge to a constant < 1)
#   Right: quad_ratio for Newton  (should converge to a constant)
# ---------------------------------------------------------------------------
warm_inits = [i for i in inits_present if i.startswith("warm")]
warm_inits = inits_present

fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(10, 4.0))
fig.subplots_adjust(left=0.09, right=0.97, bottom=0.14,
                    top=0.78, wspace=0.34)

warm_colors = [ORANGE, GREEN, GREY, RED]  # cycle through these for the warm starts
warm_markers = ["o", "s", "^"]

for i, init in enumerate(warm_inits):
    color  = warm_colors[i % len(warm_colors)]
    marker = warm_markers[i % len(warm_markers)]
    label  = INIT_LABEL.get(init, init)

    # Picard lin_ratio (skip k=0 where lin=0)
    d_p = [r for r in get(rows, "picard", init) if r["lin_ratio"] > 0]
    if d_p:
        ax0.plot([r["iter"] for r in d_p],
                 [r["lin_ratio"] for r in d_p],
                 color=color, ls="--", lw=LW,
                 marker=marker, ms=MS,
                 markerfacecolor="white", markeredgewidth=2.0,
                 label=label)

    # Newton quad_ratio (skip k=0, and drop spurious blow-ups).
    # Once |F^k| reaches machine precision (~1e-13), |F^{k+1}|/|F^k|^2
    # divides noise by a tiny squared number and explodes — these values
    # are floating-point artifacts, not the true quadratic constant.
    # Keep only iterations where the previous residual is safely above the
    # noise floor and the ratio is within a physically plausible range.
    QUAD_MAX  = 10.0     # plausible upper bound for the quadratic constant
    RES_FLOOR = 1e-11    # below this, |F| is dominated by round-off
    d_n = [r for r in get(rows, "newton", init)
           if r["quad_ratio"] > 0
           and r["quad_ratio"] < QUAD_MAX
           and r["residual_norm"] > RES_FLOOR]
    if d_n:
        ax1.plot([r["iter"] for r in d_n],
                 [r["quad_ratio"] for r in d_n],
                 color=color, ls="-", lw=LW,
                 marker=marker, ms=MS,
                 markerfacecolor="white", markeredgewidth=2.0,
                 label=label)

ax0.set_xlabel("Iteration $k$", fontsize=FS)
ax0.set_ylabel(r"$\|F^{k+1}\| / \|F^k\|$ (lin ratio)", fontsize=FS)
ax0.set_title("Picard — linear convergence rate", fontsize=FS)
ax0.tick_params(labelsize=TS)
ax0.legend(fontsize=TS, framealpha=0.92, loc="upper center",
           bbox_to_anchor=(0.5, 1.32), ncol=3)
ax0.grid(True, ls="--", lw=0.4, alpha=0.5)
ax0.set_ylim(bottom=0)

ax1.set_xlabel("Iteration $k$", fontsize=FS)
ax1.set_ylabel(r"$\|F^{k+1}\| / \|F^k\|^2$ (quad ratio)", fontsize=FS)
ax1.set_title("Newton — quadratic convergence rate", fontsize=FS)
ax1.tick_params(labelsize=TS)
ax1.legend(fontsize=TS, framealpha=0.92, loc="upper center",
           bbox_to_anchor=(0.5, 1.32), ncol=3)
ax1.grid(True, ls="--", lw=0.4, alpha=0.5)
ax1.set_ylim(bottom=0)

fig.savefig("exp10_plot_B.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved exp10_plot_B.png")