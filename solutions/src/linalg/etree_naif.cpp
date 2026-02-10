#include <vector>

struct CSC {
    int n;
    std::vector<int> debut_colonne;  
    std::vector<int> ligne_indice;  
};

std::vector<int> elimination_tree(const CSC& A)
{
    int n = A.n;
    std::vector<int> parent(n, -1);   // on ne compte pas le dernier sommet, k < j donc k varie dans 1..n-1
    std::vector<int> ancetre(n, -1);
    for (int k = 0; k < n; ++k) {
        for(int id = A.debut_colonne[k]; id < A.debut_colonne[k + 1]; ++id){ // On selectionne seulement les élements non nuls de la colonne courante
            int j = A.ligne_indice[id];
            while (j != -1 && j != k) {
                    int next = ancetre[j];
                    ancetre[j] = k;

                    if (next == -1) {
                        parent[j] = k;
                        break;
                    }
                    j = next;
                }
            }
        }
    return parent;

}