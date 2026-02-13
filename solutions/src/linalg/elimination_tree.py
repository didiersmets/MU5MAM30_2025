import numpy as np
from scipy.sparse import csr_matrix

def elimination_tree(A_csr, assume_symmetric=True,show_ancestor = False):
    """
    Fonction qui prend en entrée une matrice A symétrique et qui renvoie un vecteur de parent 
    issu de l'elimination tree.
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
    if show_ancestor == True:
        return parent, ancestor
    return parent

# test :
A_dense = np.array([
    [4, 1, 0, 0],
    [1, 5, 2, 3],
    [0, 2, 6, 0],
    [0, 3, 0, 7]
], dtype=float)

A = csr_matrix(A_dense)
print(elimination_tree(A))

A_dense2 = np.array([
    [6, 1, 0, 0, 0, 0],
    [1, 7, 2, 3, 0, 0],
    [0, 2, 5, 0, 0, 0],
    [0, 3, 0, 6, 4, 0],
    [0, 0, 0, 4, 8, 2],
    [0, 0, 0, 0, 2, 5]
], dtype=float)

A2 = csr_matrix(A_dense2)
print(elimination_tree(A2))