
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags, kron, eye, csr_matrix
from cholesky import sparse_cholesky_up_looking, elimination_tree, construction_S_L

def dense_spd(n, eps=1e-1, seed=0):
    rng = np.random.default_rng(seed)
    M = rng.standard_normal((n, n))
    A = M @ M.T + eps * np.eye(n)
    return A

def Laplacien(n):
     return diags(
        [-np.ones(n-1), 2*np.ones(n), -np.ones(n-1)],
        [-1, 0, 1],
        format="csr"
    )     

from scipy.sparse import diags, kron, eye

def Laplacien_2D(m):
    """
    Laplacien 2D sur une grille m x m.
    Taille finale : n = m^2
    """
    # Laplacien 1D
    T = diags(
        [-np.ones(m-1), 2*np.ones(m), -np.ones(m-1)],
        [-1, 0, 1],
        format="csr"
    )

    I = eye(m, format="csr")

    # Produit de Kronecker
    A = kron(I, T) + kron(T, I)

    return A


def s_per_row_from_L(L):
    # L est supposée dense (np.array) ou sparse CSR; on renvoie s_i pour i=0..n-1
    n = L.shape[0]
    s = np.zeros(n, dtype=int)

    # version dense (simple)
    if isinstance(L, np.ndarray):
        for i in range(n):
            s[i] = np.count_nonzero(L[i, :i])  # strictement sous-diagonale
        return s

    # version sparse CSR
    L = L.tocsr()
    for i in range(n):
        row_idx = L.indices[L.indptr[i]:L.indptr[i+1]]
        s[i] = np.sum(row_idx < i)
    return s

def ratio_theorique_dense(ns, ops):
        ns = np.array(ns, dtype=float)
        ops = np.array(ops, dtype=float)
        return ops / (ns**3 / 3.0)

def ratio_theorique_sparse(sum_si2, ops):
    sum_si2 = np.array(sum_si2, dtype=float)
    ops = np.array(ops, dtype=float)
    return ops / sum_si2


if __name__ == "__main__":

    print("TEST 1 : Matrice 4x4 , exemple du rapport : ")
    print("Test de l'elimination tree ")
    A_dense_1 = np.array([
        [4, 1, 0, 0],
        [1, 5, 2, 3],
        [0, 2, 6, 0],
        [0, 3, 0, 7]
    ], dtype=float)

    A1 = csr_matrix(A_dense_1)
    Parents_A1 = elimination_tree(A1,True)
    print("Parents associés à A :" ,Parents_A1)
    print("Liste des éléments non nuls de L :",construction_S_L(3,A1,Parents_A1))

    print(" " * 70)
    print("Test de la factorisation up-looking Cholesky")

    L_dense = sparse_cholesky_up_looking(A_dense_1,True)[0]
    nb_operations = sparse_cholesky_up_looking(A_dense_1,True)[1]
    print("L =")
    print(L_dense)
    print("Nombres d'opérations :", nb_operations)
    print("Vérification A ≈ L L^T :")
    print("||A - L L^T|| =", np.linalg.norm(A_dense_1 - L_dense @ L_dense.T))

    sum_si2_list = []
    ops_s = []
    #ns_s = [50,100,200,300,500,700]  
    ns_s = [20,30,50,70,100]  
    for n in ns_s:
        A = Laplacien_2D(n)
        L, ops = sparse_cholesky_up_looking(A, return_opcount=True)
        s = s_per_row_from_L(L)
        sum_si2_list.append(np.sum(s**2))
        ops_s.append(ops)
    
    plt.figure()
    plt.loglog(sum_si2_list,ops_s,  marker="o", label="Sparse")
    plt.xlabel(r"s_i^2")
    plt.ylabel("Nombre d'opérations")
    plt.grid(True)
    plt.legend()
    plt.savefig("ops_vs_n_sparse.pdf")
    plt.close()


    plt.figure()
    plt.semilogx(ns_s, ratio_theorique_sparse(sum_si2_list, ops_s), marker="o")
    plt.xlabel("n")
    plt.ylabel(r"ops / $\sum_i s_i^2$")
    plt.grid(True)
    plt.savefig("ratio_vs_n_sparse.pdf")
    plt.close()

    plt.figure()
    plt.semilogx(sum_si2_list, ratio_theorique_sparse(sum_si2_list, ops_s), marker="o")
    plt.xlabel("sum")
    plt.ylabel(r"ops / $\sum_i s_i^2$")
    plt.grid(True)
    plt.savefig("ratio_vs_s_i2_sparse.pdf")
    plt.close()

    print("fin")
    print(sum_si2_list)
    print(ns_s)

    # --- Dense ---
    # ns_d = [50,100,200,300,500,700]   
    # ops_d = []
    # for n in ns_d:
    #     A = dense_spd(n)
    #     _, ops = sparse_cholesky_up_looking(A, return_opcount=True)
    #     ops_d.append(ops)

    # plt.figure()
    # plt.loglog(ns_d,  ops_d,  marker="o", label="SPD dense")
    # plt.xlabel("n")
    # plt.ylabel("Nombre d'opérations")
    # plt.grid(True)
    # plt.legend()
    # plt.savefig("ops_vs_n.pdf")
    # plt.close()

    # # Figure 3 : ratio ops / (n^3/3)

    # plt.figure()
    # plt.loglog(ns_d,  ratio(ns_d,  ops_d),  marker="o", label="Dense")
    # plt.xlabel("n")
    # plt.ylabel(r"ops / (n^3/3)")
    # plt.grid(True)
    # plt.legend()
    # plt.savefig("ratio_vs_n.pdf")
    # plt.close()
    # print("fin")

   