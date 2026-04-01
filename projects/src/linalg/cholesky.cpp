#include "../../include/linalg/cholesky.h"
#include <vector>
#include <cassert>

// Algorithme d'elimination tree (path compression)
void elimination_tree(const CSRMatrix& A, std::vector<int>& parent) {
    size_t n = A.rows;
    parent.assign(n, -1); // -1 = racine
    std::vector<int> ancestor(n, -1);

    for (size_t i = 0; i < n; ++i) {
        ancestor[i] = -1;
        // Parcours des colonnes j < i dans la ligne i (sous-diagonale)
        for (size_t k = A.row_start[i]; k < A.row_start[i+1]; ++k) {
            int j = A.col[k];
            if (j >= (int)i) continue;
            int jroot = j;
            // Path compression
            while (ancestor[jroot] != -1 && ancestor[jroot] != (int)i) {
                int l = ancestor[jroot];
                ancestor[jroot] = (int)i;
                jroot = l;
            }
            if (ancestor[jroot] == -1) {
                ancestor[jroot] = (int)i;
                parent[jroot] = (int)i;
            }
        }
    }
}
