import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse import issparse

# Algorithme Elimination Tree (Pseudo-code 4.2 du livre)

def elimination_tree(A_csr,show_ancestor=False):
    """
    Fonction qui prend en entrée une matrice A symétrique et qui renvoie un vecteur de parent 
    issu de l'Elimination Tree
    """
    n = A_csr.shape[0]
    assert A_csr.shape[0] == A_csr.shape[1] # A doit être carré 

    # A l'initialisation, il n'y a pas de parent ni d'ancêtre. 
    # Comme -1 est un indice non valide par rapport à 0, on initialise à -1 chaque coefficient
    # Précision en float64 pour de gros calculs et grosses matrices : assurer un bon rapport coût/précision 
    parent = -np.ones(n, dtype=np.int64) 
    ancestor = -np.ones(n, dtype=np.int64)

    for i in range(n):
        ancestor[i] = -1  
        row_start = A_csr.indptr[i] # Structure propre des matrices CSR
        row_end = A_csr.indptr[i + 1]
        neigh = A_csr.indices[row_start:row_end] # Voisins de i dans le graphe de A
        for j in neigh:
            if j >= i :
                continue
            jroot = j
            while ancestor[jroot] != -1 and ancestor[jroot] != i:
                l = ancestor[jroot]
                ancestor[jroot] = i
                jroot = l
            if ancestor[jroot] == -1:
                ancestor[jroot] = i
                parent[jroot] = i
    if show_ancestor:
        print("Ancêtres associés à A : ", ancestor)
    return parent

# Phase symbolique locale : détermination de S{L_{i, 1:i-1}} (pseudo-code 4.3 du livre)

def construction_S_L(i, A_csr, parent):
    """
    Fonction qui renvoie la liste triée des colonnes k<i où L[i,k] est non nul
    (S(L_{i,1:i-1}) dans le pseudo-code)
    """
    start = A_csr.indptr[i]
    end = A_csr.indptr[i+1]
    neigh = A_csr.indices[start:end]

    L = set() # Ensemble des colonnes déjà trouvées pour la ligne i de L
    for j in neigh: # On regarde tous les voisins j de i dans le graphe de A
        if j >= i:  # On ne garde que les colonnes plus petites que i (car L est triangulaire inférieur)
            continue
        v = int(j)
        while v != -1 and v < i and v not in L:
            L.add(v) # On ajoute v dans la structure non nulle de la ligne i
            v = int(parent[v])
    # On renvoie la liste triée des colonnes non nulles de L[i, :]
    return sorted(L)

def A_coefficient(A, i, j):
    """
    Fonction qui renvoie le coefficient A[i,j] d'une matrice A, pour une matrice sparse/csr ou dense
    """
    if issparse(A):
        # On définit la structure csr
        row_start = A.indptr[i]
        row_end = A.indptr[i + 1]
        vals  = A.data[row_start:row_end]
        cols = A.indices[row_start:row_end]

        # Vérifier si la colonne j fait bien partie des colonnes non nulles de la ligne i
        idx = np.searchsorted(cols, j)
        if idx < len(cols) and cols[idx] == j:
            return vals[idx]
        else:
            return 0.0
    else:
        return A[i, j]


def sparse_cholesky_up_looking(A,return_opcount=False):
    """
    Factorisation de Cholesky creuse up-looking A = LL^T (pseudo-code 5.7 du livre)
    """

    if issparse(A):
        n = A.shape[0]
        assert A.shape[0] == A.shape[1], "A doit être carrée"
    else:
        n = A.shape[0]
        assert A.shape[0] == A.shape[1], "A doit être carrée"

    # Compteur d'opérations
    op_count = 0

    # On convertit toujours en CSR pour le calcul de l'elimination tree car notre fonction requiert le format CSR 
    if issparse(A):
        A_csr = A
    else:
        A_csr = csr_matrix(A)

    parent = elimination_tree(A_csr)

    # Initialisation de L format dense : 
    L = np.zeros((n, n), dtype=np.float64)
    a_00 = A_coefficient(A, 0, 0)
    assert a_00 > 0, "La matrice n'est pas définie positive (a_{00} <= 0)"
    L[0, 0] = np.sqrt(a_00)

    for i in range(1, n):
        # Phase symbolique 
        S = construction_S_L(i, A_csr, parent)
        # Initialisation du vecteur solution du système linéaire Ly = b 
        y = np.zeros(i, dtype=np.float64)
        for k in S:
            # b[k] = A[i, k] 
            b_k = A_coefficient(A, k, i)
            s = 0.0
            for m in S:
                if m >= k:
                    break  # S est trié, donc on s'arrête dès que m >= k
                s += L[k, m] * y[m]
                op_count += 2 # 1 multiplication + 1 addition

            assert L[k, k] != 0.0, f"Pivot nul en position ({k},{k}), la matrice n'est pas SPD"
            y[k] = (b_k - s) / L[k, k]
            op_count += 2  # 1 soustraction + 1 division

        for k in S:
            L[i, k] = y[k]


        a_ii = A_coefficient(A, i, i)
        s = 0.0
        for k in S:
            s += L[i, k] ** 2
            op_count += 2  # 1 multiplication + 1 addition

        l_ii = a_ii - s
        op_count += 1  # 1 soustraction
        L[i, i] = np.sqrt(l_ii)
    if return_opcount:
        return L, op_count

    return L