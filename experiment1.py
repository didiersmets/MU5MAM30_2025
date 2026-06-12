"""
plot_experiment_1.py  —  Experiment 1: Convergence analysis

Two PNGs:
  exp1_plots_scherk.png  — Scherk (test case 3): errors vs h + convergence rates
  exp1_plots_degenerate.png — Linear (tc1) and flat (tc2): errors vs h

Usage:  python plot_experiment_1.py [--csv PATH]
"""

import argparse, csv
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--csv", default="experiment_1_convergence.csv")
args = parser.parse_args()

# ---------------------------------------------------------------------------
# Load all test cases
# ---------------------------------------------------------------------------
all_rows = []
with open(args.csv, newline="") as fh:
    for r in csv.DictReader(fh):
        all_rows.append({k: float(r[k]) for k in r})

def get_case(tc):
    rows = [r for r in all_rows if int(r["test_case"]) == tc]
    return sorted(rows, key=lambda r: r["h"], reverse=True)

tc1 = get_case(1)
tc2 = get_case(2)
tc3 = get_case(3)

BLUE   = "#2166ac"
RED    = "#d6604d"
GREEN  = "#4dac26"
ORANGE = "#f4a742"
GREY   = "#555555"
FS = 12
TS = 10
LW = 2.8
MS = 8

# ===========================================================================
# Figure 1: Scherk (test case 3)
# ===========================================================================
rows  = tc3
h     = np.array([r["h"]        for r in rows])
L2    = np.array([r["L2_error"] for r in rows])
H1    = np.array([r["H1_error"] for r in rows])
rL2   = np.array([r["rate_L2"]  for r in rows])
rH1   = np.array([r["rate_H1"]  for r in rows])
levels = np.arange(len(rows))

L2_ref = L2[0] * (h / h[0])**2
H1_ref = H1[0] * (h / h[0])**1

fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(9, 3.6))
fig.subplots_adjust(left=0.10, right=0.97, bottom=0.17, top=0.87, wspace=0.38)

ax0.loglog(h, L2_ref, color=BLUE, ls=":", lw=2.2, label=r"$O(h^2)$")
ax0.loglog(h, H1_ref, color=RED,  ls=":", lw=2.2, label=r"$O(h)$")
ax0.loglog(h, L2, color=BLUE, ls="--", marker="o", ms=MS, lw=LW,
           markerfacecolor="white", markeredgewidth=2.2, label=r"$L^2$ error")
ax0.loglog(h, H1, color=RED,  ls="--", marker="s", ms=MS, lw=LW,
           markerfacecolor="white", markeredgewidth=2.2, label=r"$H^1$ error")
ax0.invert_xaxis()
ax0.set_xlabel(r"Mesh size $h$", fontsize=FS)
ax0.set_ylabel("Error", fontsize=FS)
ax0.set_title(r"$L^2$ and $H^1$ errors vs $h$", fontsize=FS)
ax0.tick_params(labelsize=TS)
ax0.legend(fontsize=TS, framealpha=0.92)
ax0.grid(True, which="both", ls="--", lw=0.4, alpha=0.5)

ax1.axhline(2.0, color=BLUE, ls=":", lw=2.2, label=r"target $p=2$")
ax1.axhline(1.0, color=RED,  ls=":", lw=2.2, label=r"target $p=1$")
ax1.plot(levels[1:], rL2[1:], color=BLUE, ls="--", marker="o", ms=MS, lw=LW,
         markerfacecolor="white", markeredgewidth=2.2, label=r"rate $L^2$")
ax1.plot(levels[1:], rH1[1:], color=RED,  ls="--", marker="s", ms=MS, lw=LW,
         markerfacecolor="white", markeredgewidth=2.2, label=r"rate $H^1$")
ax1.set_xlabel("Refinement level", fontsize=FS)
ax1.set_ylabel(r"Convergence rate $\hat{p}$", fontsize=FS)
ax1.set_title("Empirical convergence rates", fontsize=FS)
ax1.set_xticks(levels[1:])
ax1.set_xticklabels([f"$N={int(r['N_subdiv'])}$" for r in rows[1:]], fontsize=TS)
ax1.tick_params(labelsize=TS)
ax1.legend(fontsize=TS, framealpha=0.92)
ax1.grid(True, ls="--", lw=0.4, alpha=0.5)
ax1.set_ylim(0.8, 2.2)

#fig.suptitle(r"Experiment 1 — Scherk minimal surface (test case 3)", fontsize=FS+1)
fig.savefig("exp1_plots_scherk.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved exp1_plots_scherk.png")

# ===========================================================================
# Figure 2: degenerate cases (test cases 1 and 2)
# Linear and flat solutions are P1-exact: all errors identically zero.
# Plot: bar chart of error values at each mesh level (all zero = confirmed).
# ===========================================================================
fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(9, 3.6))
fig.subplots_adjust(left=0.10, right=0.97, bottom=0.22, top=0.87, wspace=0.38)

# ── Left: linear u = x+y (test case 1) ──────────────────────────────────────
h1   = np.array([r["h"]         for r in tc1])
L2_1 = np.array([r["L2_error"]  for r in tc1])
H1_1 = np.array([r["H1_error"]  for r in tc1])
N1   = [int(r["N_subdiv"]) for r in tc1]
x1   = np.arange(len(tc1))

OFFSET = 0.004   # small visual offset so the two zero lines don't overlap
ax0.plot(x1, L2_1 + OFFSET, color=BLUE, ls="--", marker="o", ms=MS, lw=LW,
         markerfacecolor="white", markeredgewidth=2.2, label=r"$L^2$ error")
ax0.plot(x1, H1_1,           color=RED,  ls="--", marker="s", ms=MS, lw=LW,
         markerfacecolor="white", markeredgewidth=2.2, label=r"$H^1$ error")
ax0.set_xticks(x1)
ax0.set_xticklabels([f"$N={n}$" for n in N1], fontsize=TS)
ax0.set_xlabel("Mesh refinement", fontsize=FS)
ax0.set_ylabel("Error", fontsize=FS)
ax0.set_title(r"Test case 1: linear $u = x+y$", fontsize=FS)
ax0.tick_params(labelsize=TS)
ax0.legend(fontsize=TS, framealpha=0.92)
ax0.set_ylim(-0.02, 0.06)
ax0.grid(True, ls="--", lw=0.4, alpha=0.5)
ax0.text(0.5, 0.72, "All errors $= 0$\n(P1-exact solution)",
         transform=ax0.transAxes, ha="center", va="center",
         fontsize=TS+1, color=GREY,
         bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=GREY, alpha=0.85))

# ── Right: flat u = 0 on disk (test case 2) ─────────────────────────────────
h2   = np.array([r["h"]         for r in tc2])
L2_2 = np.array([r["L2_error"]  for r in tc2])
H1_2 = np.array([r["H1_error"]  for r in tc2])
N2   = [int(r["N_subdiv"]) for r in tc2]
x2   = np.arange(len(tc2))

ax1.plot(x2, L2_2 + OFFSET, color=BLUE, ls="--", marker="o", ms=MS, lw=LW,
         markerfacecolor="white", markeredgewidth=2.2, label=r"$L^2$ error")
ax1.plot(x2, H1_2,           color=RED,  ls="--", marker="s", ms=MS, lw=LW,
         markerfacecolor="white", markeredgewidth=2.2, label=r"$H^1$ error")
ax1.set_xticks(x2)
ax1.set_xticklabels([f"$N={n}$" for n in N2], fontsize=TS)
ax1.set_xlabel("Mesh refinement", fontsize=FS)
ax1.set_ylabel("Error", fontsize=FS)
ax1.set_title(r"Test case 2: flat $u = 0$, disk $R=1$", fontsize=FS)
ax1.tick_params(labelsize=TS)
ax1.legend(fontsize=TS, framealpha=0.92)
ax1.set_ylim(-0.02, 0.06)
ax1.grid(True, ls="--", lw=0.4, alpha=0.5)
ax1.text(0.5, 0.72, "All errors $= 0$\n(P1-exact solution)",
         transform=ax1.transAxes, ha="center", va="center",
         fontsize=TS+1, color=GREY,
         bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=GREY, alpha=0.85))

#fig.suptitle("Experiment 1 — Degenerate test cases (P1-exact solutions)",
#             fontsize=FS+1)
fig.savefig("exp1_plots_degenerate.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved exp1_plots_degenerate.png")