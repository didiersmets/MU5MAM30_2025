#include "P2.h"
#include "adjacency.h"
#include "fem_matrix.h"
#include "mass_P2.h"
#include "mesh.h"
#include "sparse_matrix.h"
#include "stiffness_P2.h"

#include <map>
#include <utility>


// assigns a global index to every edge of the triangle
static void build_edge_numbering(const Mesh &m,
                                 std::map<std::pair<uint32_t,uint32_t>, uint32_t> &edge_id,
                                 uint32_t &edge_count)
{
    edge_id.clear();
    edge_count = 0;

    const TArray<uint32_t> &idx = m.indices;
    size_t tri_count = m.triangle_count();

    for (size_t t = 0; t < tri_count; ++t) {
        uint32_t a = idx[3*t+0];
        uint32_t b = idx[3*t+1];
        uint32_t c = idx[3*t+2];

        // edges of the triangle (ordered min,max para to stay consistent)
        std::pair<uint32_t,uint32_t> e0(std::min(a,b), std::max(a,b));
        std::pair<uint32_t,uint32_t> e1(std::min(b,c), std::max(b,c));
        std::pair<uint32_t,uint32_t> e2(std::min(c,a), std::max(c,a));

        if (!edge_id.count(e0)) edge_id[e0] = edge_count++;
        if (!edge_id.count(e1)) edge_id[e1] = edge_count++;
        if (!edge_id.count(e2)) edge_id[e2] = edge_count++;
    }
}

// Builds the pattern for the mass and rigid matrix
void build_P2_CSRPattern(const Mesh &m, CSRPattern &P)
{
    size_t vtx_count = m.vertex_count();
    size_t tri_count = m.triangle_count();

    // 1. Giving an index to vertices
    std::map<std::pair<uint32_t,uint32_t>, uint32_t> edge_id;
    uint32_t edge_count = 0;
    build_edge_numbering(m, edge_id, edge_count);

    size_t ndof = vtx_count + edge_count;

    P.symmetric = true;
    P.rows = P.cols = ndof;

    // 2. List of columns for every row
    std::vector<std::vector<uint32_t>> rows(ndof);

    const TArray<uint32_t> &idx = m.indices;
    for (size_t t = 0; t < tri_count; ++t) {
        uint32_t a = idx[3*t+0];
        uint32_t b = idx[3*t+1];
        uint32_t c = idx[3*t+2];

        uint32_t dof[6];
        // vertex
        dof[0] = a;
        dof[1] = b;
        dof[2] = c;

        // aristas -> índices globales desplazados
        auto e0 = std::make_pair(std::min(a,b), std::max(a,b));
        auto e1 = std::make_pair(std::min(b,c), std::max(b,c));
        auto e2 = std::make_pair(std::min(c,a), std::max(c,a));
        dof[3] = vtx_count + edge_id[e0];
        dof[4] = vtx_count + edge_id[e1];
        dof[5] = vtx_count + edge_id[e2];

        // 3. connect every local dof  (i >= j for simetry)
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j <= i; ++j) {
                uint32_t I = dof[i];
                uint32_t J = dof[j];
                rows[I].push_back(J);
            }
        }
    }

    // 4. Eliminar duplicados y construir CSRPattern
    P.row_start.resize(ndof + 1);
    size_t nnz = 0;
    for (size_t i = 0; i < ndof; ++i) {
        auto &r = rows[i];
        std::sort(r.begin(), r.end());
        r.erase(std::unique(r.begin(), r.end()), r.end());
        P.row_start[i] = nnz;
        nnz += r.size();
    }
    P.row_start[ndof] = nnz;

    P.col.resize(nnz);
    for (size_t i = 0; i < ndof; ++i) {
        auto &r = rows[i];
        size_t offset = P.row_start[i];
        for (size_t k = 0; k < r.size(); ++k) {
            P.col[offset + k] = r[k];
        }
    }
}

void build_P2_mass_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &M)
{
    size_t vtx_count = m.vertex_count();
    size_t tri_count = m.triangle_count();

    // numerar aristas igual que en el patrón
    std::map<std::pair<uint32_t,uint32_t>, uint32_t> edge_id;
    uint32_t edge_count = 0;
    build_edge_numbering(m, edge_id, edge_count);

    size_t ndof = vtx_count + edge_count;

    M.symmetric = true;
    M.rows = M.cols = ndof;
    M.nnz = P.col.size();
    M.row_start = P.row_start.data();
    M.col = P.col.data();
    M.data.assign(M.nnz, 0.0);

    const TArray<uint32_t> &idx = m.indices;
    for (size_t t = 0; t < tri_count; ++t) {
        uint32_t a = idx[3*t+0];
        uint32_t b = idx[3*t+1];
        uint32_t c = idx[3*t+2];

        Vec3f A = m.positions[a];
        Vec3f B = m.positions[b];
        Vec3f C = m.positions[c];
        Vec3d AB = { (double)B[0] - (double)A[0],
                     (double)B[1] - (double)A[1],
                     (double)B[2] - (double)A[2] };
        Vec3d AC = { (double)C[0] - (double)A[0],
                     (double)C[1] - (double)A[1],
                     (double)C[2] - (double)A[2] };

        double Mloc[36];
        mass_P2(AB, AC, Mloc);

        uint32_t dof[6];
        dof[0] = a;
        dof[1] = b;
        dof[2] = c;
        auto e0 = std::make_pair(std::min(a,b), std::max(a,b));
        auto e1 = std::make_pair(std::min(b,c), std::max(b,c));
        auto e2 = std::make_pair(std::min(c,a), std::max(c,a));
        dof[3] = vtx_count + edge_id[e0];
        dof[4] = vtx_count + edge_id[e1];
        dof[5] = vtx_count + edge_id[e2];

        // ensamblaje: M(dof[i], dof[j]) += Mloc[i*6 + j]
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j <= i; ++j) { // simetría
                uint32_t I = dof[i];
                uint32_t J = dof[j];
                double val = Mloc[i*6 + j];
                M(I, J) += val;
            }
        }
    }
}


void build_P2_stiffness_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &S)
{
    size_t vtx_count = m.vertex_count();
    size_t tri_count = m.triangle_count();

    std::map<std::pair<uint32_t,uint32_t>, uint32_t> edge_id;
    uint32_t edge_count = 0;
    build_edge_numbering(m, edge_id, edge_count);

    size_t ndof = vtx_count + edge_count;

    S.symmetric = true;
    S.rows = S.cols = ndof;
    S.nnz = P.col.size();
    S.row_start = P.row_start.data();
    S.col = P.col.data();
    S.data.assign(S.nnz, 0.0);

    const TArray<uint32_t> &idx = m.indices;
    for (size_t t = 0; t < tri_count; ++t) {
        uint32_t a = idx[3*t+0];
        uint32_t b = idx[3*t+1];
        uint32_t c = idx[3*t+2];

        Vec3f A = m.positions[a];
        Vec3f B = m.positions[b];
        Vec3f C = m.positions[c];
        Vec3d AB = { (double)B[0] - (double)A[0],
                     (double)B[1] - (double)A[1],
                     (double)B[2] - (double)A[2] };
        Vec3d AC = { (double)C[0] - (double)A[0],
                     (double)C[1] - (double)A[1],
                     (double)C[2] - (double)A[2] };

        double Sloc[36];
        stiffness_P2(AB, AC, Sloc);

        uint32_t dof[6];
        dof[0] = a;
        dof[1] = b;
        dof[2] = c;
        auto e0 = std::make_pair(std::min(a,b), std::max(a,b));
        auto e1 = std::make_pair(std::min(b,c), std::max(b,c));
        auto e2 = std::make_pair(std::min(c,a), std::max(c,a));
        dof[3] = vtx_count + edge_id[e0];
        dof[4] = vtx_count + edge_id[e1];
        dof[5] = vtx_count + edge_id[e2];

        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j <= i; ++j) {
                uint32_t I = dof[i];
                uint32_t J = dof[j];
                double val = Sloc[i*6 + j];
                S(I, J) += val;
            }
        }
    }
}
