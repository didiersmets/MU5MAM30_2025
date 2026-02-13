import numpy as np
from elimination_tree import *

def sorted_list_non_zeros_L(i, A_csr, parent):
    """
    Fonction qui renvoie la liste triée des colonnes k<i où L[i,k] est non nul
    (S(L_{i,1:i-1}) dans le pseudo-code)
    """
    start = A_csr.indptr[i]
    end = A_csr.indptr[i+1]
    neigh = A_csr.indices[start:end]

    lst = set() # Ensemble des colonnes déjà trouvées pour la ligne i de L
    for j in neigh: # On regarde tous les voisins j de i dans le graphe de A
        if j >= i:  # On ne garde que les colonnes plus petites que i (car L est triangulaire inférieur)
            continue
        v = int(j)
        while v != -1 and v < i and v not in lst:
            lst.add(v) # On ajoute v dans la structure non nulle de la ligne i
            v = int(parent[v])
    # On renvoie la liste triée des colonnes non nulles de L[i, :]
    return sorted(lst)

print(sorted_list_non_zeros_L(3,A,elimination_tree(A)))