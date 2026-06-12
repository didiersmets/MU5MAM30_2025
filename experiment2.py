"""
plot_experiment_2.py  —  Experiment 2: Newton convergence (cold + warm starts)

Produces a single PNG with three panels side by side:
  Left   — |F^(k)| vs iteration
  Centre — linear ratio |F^{k+1}|/|F^k|
  Right  — quadratic ratio |F^{k+1}|/|F^k|^2

Usage:
    python plot_experiment_2.py [--csv PATH]
"""

import argparse, csv
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--csv", default="experiment_2_newton_convergence.csv")
args = parser.parse_args()

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
rows = []
with open(args.csv, newline="") as fh:
    for r in csv.DictReader(fh):
        entry = {"init": r["init"]}
        entry.update({k: float(r[k]) for k in
                      ("iter","residual_norm","correction_norm",
                       "lin_ratio","quad_ratio","loglog_slope")})
        rows.append(entry)

groups = defaultdict(list)
for r in rows:
    groups[r["init"]].append(r)

def arr(g, key):
    return np.array([r[key] for r in g])

ORDER = ["cold", "warm_1e-01", "warm_1e-02", "warm_1e-04"]
KEYS  = [k for k in ORDER if k in groups]

STYLES = {
    "cold":       dict(color="#444444", ls="--", marker="D",  ms=9,  lw=3.0,
                       label="Cold start (interior $= 0$)"),
    "warm_1e-01": dict(color="#d6604d", ls="-",  marker="o",  ms=9,  lw=3.0,
                       label=r"Warm $\varepsilon = 10^{-1}$"),
    "warm_1e-02": dict(color="#2166ac", ls="-",  marker="s",  ms=9,  lw=3.0,
                       label=r"Warm $\varepsilon = 10^{-2}$"),
    "warm_1e-04": dict(color="#4dac26", ls="-",  marker="^",  ms=9,  lw=3.0,
                       label=r"Warm $\varepsilon = 10^{-4}$"),
}
GREY = "#666666"
FS   = 13   # base font size for labels / titles
TS   = 11   # tick label size

# ---------------------------------------------------------------------------
# Single figure, three panels
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
fig.subplots_adjust(left=0.07, right=0.97, bottom=0.15, top=0.88, wspace=0.35)

# ── Panel 1: residual norm ───────────────────────────────────────────────────
ax = axes[0]
for key in KEYS:
    g   = groups[key]
    itr = arr(g, "iter")
    res = arr(g, "residual_norm")
    st  = STYLES[key]
    ax.semilogy(itr, res, color=st["color"], ls=st["ls"],
                marker=st["marker"], ms=st["ms"], lw=st["lw"])

ax.axhline(1e-13, color="black", lw=1.5, ls="-.", alpha=0.5)
ax.set_xlabel("Newton iteration $k$", fontsize=FS)
ax.set_ylabel(r"$\|F^{(k)}\|_2$", fontsize=FS)
ax.set_title("Residual norm", fontsize=FS)
ax.tick_params(labelsize=TS)
ax.grid(True, which="major", ls="--", lw=0.5, alpha=0.5)
ax.grid(True, which="minor", ls=":", lw=0.25, alpha=0.3)

# ── Panel 2: linear contraction ratio ───────────────────────────────────────
ax = axes[1]
for key in KEYS:
    g    = groups[key]
    itr  = arr(g, "iter")[1:]
    lin  = arr(g, "lin_ratio")[1:]
    mask = lin > 0
    st   = STYLES[key]
    ax.semilogy(itr[mask], lin[mask], color=st["color"], ls=st["ls"],
                marker=st["marker"], ms=st["ms"], lw=st["lw"])

ax.set_xlabel("Newton iteration $k$", fontsize=FS)
ax.set_ylabel(r"$\|F^{(k+1)}\| / \|F^{(k)}\|$", fontsize=FS)
ax.set_title("Linear ratio", fontsize=FS)
ax.tick_params(labelsize=TS)
ax.grid(True, which="major", ls="--", lw=0.5, alpha=0.5)
ax.grid(True, which="minor", ls=":", lw=0.25, alpha=0.3)

# ── Panel 3: quadratic ratio ─────────────────────────────────────────────────
ax = axes[2]
for key in KEYS:
    g    = groups[key]
    itr  = arr(g, "iter")[1:]
    qr   = arr(g, "quad_ratio")[1:]
    res  = arr(g, "residual_norm")[1:]
    mask = (qr > 0) & (res > 5e-15)
    st   = STYLES[key]
    ax.semilogy(itr[mask], qr[mask], color=st["color"], ls=st["ls"],
                marker=st["marker"], ms=st["ms"], lw=st["lw"])

ax.axhline(1.0, color=GREY, lw=2.0, ls="--", alpha=0.75,
           label=r"$O(1)$ reference")
ax.set_xlabel("Newton iteration $k$", fontsize=FS)
ax.set_ylabel(r"$\|F^{(k+1)}\| / \|F^{(k)}\|^2$", fontsize=FS)
ax.set_title("Quadratic ratio", fontsize=FS)
ax.tick_params(labelsize=TS)
ax.grid(True, which="major", ls="--", lw=0.5, alpha=0.5)
ax.grid(True, which="minor", ls=":", lw=0.25, alpha=0.3)

# ── Single shared legend below all panels ────────────────────────────────────
handles = [plt.Line2D([0],[0],
                      color=STYLES[k]["color"],
                      ls=STYLES[k]["ls"],
                      marker=STYLES[k]["marker"],
                      ms=STYLES[k]["ms"],
                      lw=STYLES[k]["lw"],
                      label=STYLES[k]["label"])
           for k in KEYS]
handles.append(plt.Line2D([0],[0], color=GREY, ls="--", lw=2.0,
                           label=r"$O(1)$ reference"))

fig.legend(handles=handles, loc="lower center",
           ncol=len(handles), fontsize=11,
           framealpha=0.92, bbox_to_anchor=(0.52, -0.12))

#fig.suptitle(
#    r"Experiment 2 — Newton convergence  "
#    r"(Scherk b.c., fixed Jacobian, $16{\times}16$ mesh)",
#    fontsize=13
#)

fig.savefig("exp2_plots.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved exp2_plots.png")