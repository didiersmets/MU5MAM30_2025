import numpy as np
import matplotlib.pyplot as plt

# Cargar datos: formato "N, FEM, L2_error, H1_error"
rows = []
with open("../results/convergence.txt") as f:
    next(f)  # saltar header
    for line in f:
        parts = [p.strip() for p in line.split(",")]
        N   = int(parts[0])
        FEM = parts[1]
        L2  = float(parts[2])
        H1  = float(parts[3])
        rows.append((N, FEM, L2, H1))

# Convertir a arrays
import numpy as np
rows = np.array(rows, dtype=object)

Ns  = rows[:,0].astype(float)
FEM = rows[:,1]
L2  = rows[:,2].astype(float)
H1  = rows[:,3].astype(float)


# Separar P1 y P2
mask_P1 = (FEM == "P1")
mask_P2 = (FEM == "P2")

Ns_P1 = Ns[mask_P1]
L2_P1 = L2[mask_P1]
H1_P1 = H1[mask_P1]

Ns_P2 = Ns[mask_P2]
L2_P2 = L2[mask_P2]
H1_P2 = H1[mask_P2]

# h = 1/N
h_P1 = 1.0 / Ns_P1
h_P2 = 1.0 / Ns_P2

plt.figure(figsize=(8,6))

plt.loglog(h_P1, L2_P1, 'o-', label="P1 L2")
plt.loglog(h_P1, H1_P1, 'o--', label="P1 H1")

plt.loglog(h_P2, L2_P2, 's-', label="P2 L2")
plt.loglog(h_P2, H1_P2, 's--', label="P2 H1")

plt.xlabel("h = 1/N")
plt.ylabel("Error")
plt.title("Convergence of P1 and P2 FEM")
plt.grid(True, which="both")
plt.legend()

plt.savefig("../results/convergence_plot.png", dpi=200)
plt.show()
