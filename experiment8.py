"""
plot_experiment_8.py  —  Experiment 8: Sources of ill-conditioning (Scherk)

Three compact PNGs:
  exp8_plot_A.png  — Part A: kappa and CG vs h  (1x2 panels)
  exp8_plot_B.png  — Part B: kappa and CG vs (1+M^2)  (1x2 panels)
  exp8_plot_C.png  — Part C: C_hat vs N

Usage:  python plot_experiment_8.py [--csv PATH]
"""

import argparse, csv, math
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--csv", default="experiment_8_illconditioning.csv")
args = parser.parse_args()

rows = []
with open(args.csv, newline="") as fh:
    for r in csv.DictReader(fh):
        rows.append({k: float(v) for k, v in r.items()})

OVERFLOW  = 1e20
partA     = [r for r in rows if r["part"] == 1]
partB     = [r for r in rows if r["part"] == 2 and r["kappa"] < OVERFLOW]
partC_all = [r for r in rows if r["part"] == 3 and r["kappa"] < OVERFLOW]

BLUE   = "#2166ac"
RED    = "#d6604d"
ORANGE = "#f4a742"
GREEN  = "#4dac26"
GREY   = "#555555"
BROWN  = "#8b4513"

FS = 12    # axis label / title font size
TS = 10    # tick label size
LW = 2.8   # line width
MS = 8     # marker size

# ---------------------------------------------------------------------------
# Part A
# ---------------------------------------------------------------------------
h     = np.array([r["h"]          for r in partA])
kappa = np.array([r["kappa"]      for r in partA])
sqrtk = np.array([r["sqrt_kappa"] for r in partA])
cg    = np.array([r["cg_iters"]   for r in partA])
# anchor reference lines at coarsest point, span the full range
k_ref  = kappa[0] * (h / h[0])**(-2)
cg_ref = cg[0]    * (h / h[0])**(-1)

fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(9, 3.6))
fig.subplots_adjust(left=0.10, right=0.97, bottom=0.17, top=0.87, wspace=0.38)

ax0.loglog(h, k_ref, color=GREY, ls="--", lw=2.0,
           label=r"$\propto h^{-2}$")
ax0.loglog(h, kappa, color=BLUE, ls="-", marker="o", ms=MS, lw=LW,
           markerfacecolor="white", markeredgewidth=2.2,
           label=r"$\kappa_2(J)$")
ax0.invert_xaxis()
for i in range(1, len(partA)):
    h_mid = math.sqrt(h[i-1]*h[i])
    k_mid = math.sqrt(kappa[i-1]*kappa[i])
    slope = math.log(kappa[i]/kappa[i-1]) / math.log(h[i]/h[i-1])
    ax0.annotate(f"{slope:.2f}", xy=(h_mid, k_mid),
                 xytext=(5, 4), textcoords="offset points", fontsize=TS, color=BLUE)
ax0.set_xlabel(r"Mesh size $h$", fontsize=FS)
ax0.set_ylabel(r"$\kappa_2(J)$", fontsize=FS)
ax0.set_title(r"$\kappa_2 \sim h^{-2}$", fontsize=FS)
ax0.tick_params(labelsize=TS)
ax0.legend(fontsize=TS, framealpha=0.92)
ax0.grid(True, which="both", ls="--", lw=0.4, alpha=0.5)

ax1.loglog(h, cg_ref, color=GREY, ls="--", lw=2.0,
           label=r"$\propto h^{-1}$")
ax1.loglog(h, cg,    color=RED,  ls="-",  marker="s", ms=MS,   lw=LW,
           markerfacecolor="white", markeredgewidth=2.2, label="CG iters")
ax1.loglog(h, sqrtk, color=BLUE, ls=":",  marker="^", ms=MS-2, lw=2.0,
           label=r"$\sqrt{\kappa_2}$")
ax1.invert_xaxis()
ax1.set_xlabel(r"Mesh size $h$", fontsize=FS)
ax1.set_ylabel("CG iterations", fontsize=FS)
ax1.set_title(r"CG $\sim h^{-1}$", fontsize=FS)
ax1.tick_params(labelsize=TS)
ax1.legend(fontsize=TS, framealpha=0.92)
ax1.grid(True, which="both", ls="--", lw=0.4, alpha=0.5)

#fig.suptitle(r"Part A — Spatial penalty  [$\alpha=\pi/4$, $M=1$, vary $N$]",
#             fontsize=FS+1)
fig.savefig("exp8_plot_A.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved exp8_plot_A.png")

# ---------------------------------------------------------------------------
# Part B
# ---------------------------------------------------------------------------
M_arr  = np.array([r["M"]          for r in partB])
kappaB = np.array([r["kappa"]      for r in partB])
sqrtkB = np.array([r["sqrt_kappa"] for r in partB])
cgB    = np.array([r["cg_iters"]   for r in partB])
X      = 1 + M_arr**2
# anchor at first point, span full range
k_ref3   = kappaB[0] * (X / X[0])**(1.5)
cg_ref34 = cgB[0]    * (X / X[0])**(0.75)

fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(9, 3.6))
fig.subplots_adjust(left=0.10, right=0.97, bottom=0.17, top=0.87, wspace=0.38)

ax0.loglog(X, k_ref3, color=GREY, ls="--", lw=2.0,
           label=r"$(1{+}M^2)^{3/2}$ (theory)")
ax0.loglog(X, kappaB, color=BLUE, ls="-",  marker="o", ms=MS, lw=LW,
           markerfacecolor="white", markeredgewidth=2.2,
           label=r"$\kappa_2(J)$")
for i in range(1, len(partB)):
    x_mid = math.sqrt(X[i-1]*X[i])
    k_mid = math.sqrt(kappaB[i-1]*kappaB[i])
    slope = math.log(kappaB[i]/kappaB[i-1]) / math.log(X[i]/X[i-1])
    ax0.annotate(f"{slope:.2f}", xy=(x_mid, k_mid),
                 xytext=(5, 4), textcoords="offset points", fontsize=TS, color=BLUE)
ax0.set_xlabel(r"$1 + M^2$", fontsize=FS)
ax0.set_ylabel(r"$\kappa_2(J)$", fontsize=FS)
ax0.set_title(r"$\kappa_2$ vs $1+M^2$", fontsize=FS)
ax0.tick_params(labelsize=TS)
ax0.legend(fontsize=TS, framealpha=0.92)
ax0.grid(True, which="both", ls="--", lw=0.4, alpha=0.5)

ax1.loglog(X, cg_ref34, color=GREY, ls="--", lw=2.0,
           label=r"$(1{+}M^2)^{3/4}$ (theory)")
ax1.loglog(X, cgB,    color=RED,  ls="-",  marker="s", ms=MS,   lw=LW,
           markerfacecolor="white", markeredgewidth=2.2, label="CG iters")
ax1.loglog(X, sqrtkB, color=BLUE, ls=":",  marker="^", ms=MS-2, lw=2.0,
           label=r"$\sqrt{\kappa_2}$")
ax1.set_xlabel(r"$1 + M^2$", fontsize=FS)
ax1.set_ylabel("CG iterations", fontsize=FS)
ax1.set_title(r"CG vs $1+M^2$", fontsize=FS)
ax1.tick_params(labelsize=TS)
ax1.legend(fontsize=TS, framealpha=0.92)
ax1.grid(True, which="both", ls="--", lw=0.4, alpha=0.5)

#fig.suptitle(r"Part B — Nonlinear penalty  [$N=32$, vary $\alpha \to M$]"
#             "\n(overflow at $\\alpha=\\pi/8$ excluded)", fontsize=FS+1)
fig.savefig("exp8_plot_B.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved exp8_plot_B.png")

# ---------------------------------------------------------------------------
# Part C
# ---------------------------------------------------------------------------
by_alpha = defaultdict(list)
for r in partC_all:
    by_alpha[round(r["alpha"], 4)].append(r)

alpha_styles = {
    0.3927: (r"$\alpha=\pi/8$,   $M=0.41$",  GREEN),
    0.7854: (r"$\alpha=\pi/4$,   $M=1.00$",  BLUE),
    1.1781: (r"$\alpha=3\pi/8$,  $M=2.41$",  RED),
    1.4137: (r"$\alpha=9\pi/20$, $M=6.31$",  BROWN),
}

fig, ax = plt.subplots(figsize=(5.5, 4.0))
fig.subplots_adjust(left=0.15, right=0.97, bottom=0.14, top=0.75)

handles = []
for alpha_key in sorted(by_alpha.keys()):
    rs   = sorted(by_alpha[alpha_key], key=lambda r: r["N"])
    Ns   = np.array([r["N"] for r in rs])
    chat = np.array([r["kappa"] * r["h"]**2 / (1 + r["M"]**2)**1.5 for r in rs])
    label, color = alpha_styles.get(alpha_key, (f"α={alpha_key:.3f}", GREY))
    line, = ax.plot(Ns, chat, color=color, ls="-", marker="o",
                    ms=MS, lw=LW, label=label)
    handles.append(line)

ax.set_xlabel(r"Subdivisions $N$", fontsize=FS)
ax.set_ylabel(r"$\hat{C} = \kappa_2 h^2 / (1+M^2)^{3/2}$", fontsize=FS)
#ax.set_title(r"Part C — Combined scaling constant $\hat{C}$", fontsize=FS)
ax.set_xscale("log")
ax.set_xticks([8, 16, 32])
ax.set_xticks([], minor=True)
ax.set_xticklabels(["8", "16", "32"])
ax.set_xlim(6, 40)
ax.tick_params(labelsize=TS)
ax.grid(True, ls="--", lw=0.5, alpha=0.5)
ax.legend(fontsize=TS, framealpha=0.92, loc="upper center",
          bbox_to_anchor=(0.5, 1.32), ncol=2)

fig.savefig("exp8_plot_C.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved exp8_plot_C.png")