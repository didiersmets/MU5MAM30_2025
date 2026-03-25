#include <vector>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "P2.h"
#include "adjacency.h"
#include "adjacency_P2.h"
#include "fem_matrix.h"
#include "mass_P2.h"
#include "mesh.h"
#include "sparse_matrix.h"
#include "stiffness_P2.h"

static bool find(uint32_t x, const uint32_t *start, size_t count)
{
    for (size_t i = 0; i < count; ++i) {
       if (start[i] == x)
          return true;
    }
    return false;
}

void build_P2_CSRPattern(const Mesh &m, CSRPattern &P, const EdgeAdjacency &edges)
{
    size_t vtx_count = m.vertex_count();
    size_t tri_count = m.triangle_count();
    size_t ddl_count = vtx_count + edges.num_edges;

    P.symmetric = true;
    P.rows = P.cols = ddl_count;
    P.row_start.resize(ddl_count + 1);
    //comme P1 mais cette fois ci on prévoit 15 connexions maximum comme bcp plus de DDL
    size_t max_nnz = 15 * tri_count + ddl_count;
    P.col.resize(max_nnz);
    // Tableau temporaire pour stocker les voisins < a, initialement chaque case est vide et on va
    //remplir au fur et à mesure puis les ajouter dans P au format CSR
    std::vector<std::vector<uint32_t>> lower_adj(ddl_count);

    // on remplit l'adjacence directement pendant la création de la matrice cette fois
    for (size_t t = 0; t < tri_count; ++t) {
        uint32_t a = m.indices[3 * t + 0];
        uint32_t b = m.indices[3 * t + 1];
        uint32_t c = m.indices[3 * t + 2];
//on récupère les 6 DDL du triangle
        uint32_t ddls[6];
        ddls[0] = a;
        ddls[1] = b;
        ddls[2] = c;
        ddls[3] = vtx_count + edges.get_edge_id(a, b);
        ddls[4] = vtx_count + edges.get_edge_id(a, c);
        ddls[5] = vtx_count + edges.get_edge_id(b, c);
//Si on a pas encore ajouté la connexion on l'ajoute
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                uint32_t row = ddls[i];
                uint32_t col = ddls[j];
                // On ne garde que la partie inférieure, si on n'a pas ajouté la connexion on la rajoute
                    if (col < row && !find(col, lower_adj[row].data(), lower_adj[row].size())) {
                        lower_adj[row].push_back(col);
                    }
            }
        }
    }

    // On a toutes les connexions indicées par sommet on les met maintenant au format CSR
    size_t nnz = 0;
    for (size_t a = 0; a < ddl_count; ++a) {
        P.row_start[a] = nnz;
        for (size_t k = 0; k < lower_adj[a].size(); ++k) {
            P.col[nnz++] = lower_adj[a][k];
        }
        P.col[nnz++] = a; // La diagonale en dernière
    }
    P.row_start[ddl_count] = nnz;
    P.col.resize(nnz);
    P.col.shrink_to_fit();

    //Enfin on fait le tri par insertion exactement comme P1
    for (size_t a = 0; a < ddl_count; ++a) {
       uint32_t *__restrict to_sort = &P.col[P.row_start[a]];
       size_t count = P.row_start[a + 1] - P.row_start[a];

       for (size_t k = 1; k < count; ++k) {
          size_t j = k;
          while (j && (to_sort[j - 1] > to_sort[j])) {
             uint32_t tmp = to_sort[j - 1];
             to_sort[j - 1] = to_sort[j];
             to_sort[j] = tmp;
             j--;
          }
       }
    }
}

void build_P2_mass_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &M, const EdgeAdjacency &edges)
{
    size_t vtx_count = m.vertex_count();
    size_t tri_count = m.triangle_count();
    size_t ddl_count = P.rows; //permet de compter tout les DDL

    M.symmetric = true;
    M.rows = M.cols = ddl_count;
    M.nnz = P.col.size;
    M.row_start = P.row_start.data;
    M.col = P.col.data;
    M.data.resize(M.nnz);
    for (size_t i = 0; i < M.nnz; ++i) {
       M.data[i] = 0.0;
    }

    const TArray<uint32_t> &idx = m.indices;

    for (size_t t = 0; t < tri_count; ++t) {
       uint32_t a = idx[3 * t + 0];
       uint32_t b = idx[3 * t + 1];
       uint32_t c = idx[3 * t + 2];
       Vec3f A = m.positions[a];
       Vec3f B = m.positions[b];
       Vec3f C = m.positions[c];
       Vec3d AB = { (double)B[0] - (double)A[0],
               (double)B[1] - (double)A[1],
               (double)B[2] - (double)A[2] };
       Vec3d AC = { (double)C[0] - (double)A[0],
               (double)C[1] - (double)A[1],
               (double)C[2] - (double)A[2] };

       //On récupère les 6 indices globaux du triangle
       uint32_t dofs[6];
       dofs[0] = a;
       dofs[1] = b;
       dofs[2] = c;
       dofs[3] = vtx_count + edges.get_edge_id(a, b);
       dofs[4] = vtx_count + edges.get_edge_id(a, c);
       dofs[5] = vtx_count + edges.get_edge_id(b, c);

       //Matrice locale 6x6 au lieu de Mloc[2] avant
       double Mloc[6][6];
       mass_P2(AB, AC, Mloc);

       // même logique que P1 mais boucle car plus d'éléments
       for (int i = 0; i < 6; ++i) {
           // On s'arrête à j <= i pour ne parcourir que la partie triangulaire inférieure de Mloc
           for (int j = 0; j <= i; ++j) {

               uint32_t row = dofs[i] > dofs[j] ? dofs[i] : dofs[j];
               uint32_t col = dofs[i] > dofs[j] ? dofs[j] : dofs[i];
               M(row, col) += Mloc[i][j];
           }
       }
    }
}
void build_P2_stiffness_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &S, const EdgeAdjacency &edges)
{
    size_t vtx_count = m.vertex_count();
    size_t tri_count = m.triangle_count();
    size_t ddl_count = P.rows;

    S.symmetric = true;
    S.rows = S.cols = ddl_count;
    S.nnz = P.col.size;
    S.row_start = P.row_start.data;
    S.col = P.col.data;
    S.data.resize(S.nnz);
    for (size_t i = 0; i < S.nnz; ++i) {
       S.data[i] = 0.0;
    }

    const TArray<uint32_t> &idx = m.indices;

    for (size_t t = 0; t < tri_count; ++t) {
       // Géométrie du triangle
       uint32_t a = idx[3 * t + 0];
       uint32_t b = idx[3 * t + 1];
       uint32_t c = idx[3 * t + 2];
       Vec3f A = m.positions[a];
       Vec3f B = m.positions[b];
       Vec3f C = m.positions[c];
       Vec3d AB = { (double)B[0] - (double)A[0],
               (double)B[1] - (double)A[1],
               (double)B[2] - (double)A[2] };
       Vec3d AC = { (double)C[0] - (double)A[0],
               (double)C[1] - (double)A[1],
               (double)C[2] - (double)A[2] };

       // Les 6 indices globaux du triangle
       uint32_t dofs[6];
       dofs[0] = a;
       dofs[1] = b;
       dofs[2] = c;
       dofs[3] = vtx_count + edges.get_edge_id(a, b);
       dofs[4] = vtx_count + edges.get_edge_id(a, c);
       dofs[5] = vtx_count + edges.get_edge_id(b, c);
//Matrice 6x6 également
       double Sloc[6][6];
       stiffness_P2(AB, AC, Sloc);

       for (int i = 0; i < 6; ++i) {
           for (int j = 0; j <= i; ++j) {
               uint32_t row = dofs[i] > dofs[j] ? dofs[i] : dofs[j];
               uint32_t col = dofs[i] > dofs[j] ? dofs[j] : dofs[i];

               S(row, col) += Sloc[i][j];
           }
       }
    }
}