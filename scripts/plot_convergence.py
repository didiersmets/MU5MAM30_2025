import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# 1. Cargar datos
# ============================================================

rows = []
with open("../results/convergence.txt") as f:
    next(f)
    for line in f:
        parts = [p.strip() for p in line.split(",")]
        rows.append((int(parts[0]), parts[1], float(parts[2]), float(parts[3])))

rows = np.array(rows, dtype=object)

Ns  = rows[:,0].astype(float)
FEM = rows[:,1]
L2  = rows[:,2].astype(float)
H1  = rows[:,3].astype(float)

mask_P1 = (FEM == "P1")
mask_P2 = (FEM == "P2")

Ns_P1, L2_P1, H1_P1 = Ns[mask_P1], L2[mask_P1], H1[mask_P1]
Ns_P2, L2_P2, H1_P2 = Ns[mask_P2], L2[mask_P2], H1[mask_P2]

h_P1 = 1.0 / Ns_P1
h_P2 = 1.0 / Ns_P2

# ============================================================
# 2. Cálculo de órdenes
# ============================================================

def compute_order(h, err):
    m, b = np.polyfit(np.log(h), np.log(err), 1)
    return m

order_P1_L2 = compute_order(h_P1, L2_P1)
order_P1_H1 = compute_order(h_P1, H1_P1)
order_P2_L2 = compute_order(h_P2, L2_P2)
order_P2_H1 = compute_order(h_P2, H1_P2)

print("Orden P1 L2 =", order_P1_L2)
print("Orden P1 H1 =", order_P1_H1)
print("Orden P2 L2 =", order_P2_L2)
print("Orden P2 H1 =", order_P2_H1)

# ============================================================
# 3. Estilo académico
# ============================================================

plt.rcParams.update({
    "font.size": 13,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "figure.figsize": (7,5),
    "lines.linewidth": 1.8,
    "lines.markersize": 6,
})

fig, ax = plt.subplots()

ax.loglog(h_P1, L2_P1, 'o-', color="#1f77b4", label=f"P1 L2 (p={order_P1_L2:.2f})")
ax.loglog(h_P1, H1_P1, 'o--', color="#ff7f0e", label=f"P1 H1 (p={order_P1_H1:.2f})")

ax.loglog(h_P2, L2_P2, 's-', color="#2ca02c", label=f"P2 L2 (p={order_P2_L2:.2f})")
ax.loglog(h_P2, H1_P2, 's--', color="#d62728", label=f"P2 H1 (p={order_P2_H1:.2f})")

ax.set_xlabel("h = 1/N")
ax.set_ylabel("Error")

# ax.set_title("Convergencia de P1 y P2 FEM")

ax.legend(frameon=True, facecolor="white", edgecolor="black")

plt.tight_layout()
plt.savefig("../results/convergence_plot.pdf", dpi=300)
plt.savefig("../results/convergence_plot.png", dpi=300)
plt.show()
