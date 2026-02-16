#include "P2.h"
#include "adjacency.h"
#include "fem_matrix.h"
#include "mass_P2.h"
#include "mesh.h"
#include "sparse_matrix.h"
#include "stiffness_P2.h"

#include <set>
#include <map>
#include <utility>
#include <algorithm> // std::sort, std::unique, std::min, std::max



// CSRMatrix variants
void build_P2_CSRPattern(const Mesh &m, CSRPattern &P)
{
    size_t vtx_count = m.vertex_count();
    size_t tri_count = m.triangle_count();
    const TArray<uint32_t> &idx = m.indices;

    std::map<std::pair<uint32_t,uint32_t>, uint32_t> edge_id;
    uint32_t edge_count = 0;
    build_edge_numbering(m, edge_id, edge_count);

    size_t ndof = vtx_count + edge_count;
    P.symmetric = true;
    P.rows = P.cols = ndof;

    std::vector<std::set<uint32_t>> adj(ndof);

    for (size_t t = 0; t < tri_count; ++t) {
        uint32_t a = idx[3*t+0];
        uint32_t b = idx[3*t+1];
        uint32_t c = idx[3*t+2];

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

        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j) {
                adj[dof[i]].insert(dof[j]);
                adj[dof[j]].insert(dof[i]);
            }
    }

    P.row_start.resize(ndof + 1);
    size_t nnz = 0;
    std::vector<std::vector<uint32_t>> rows(ndof);

    for (size_t i = 0; i < ndof; ++i) {
        auto &cols = rows[i];
        cols.assign(adj[i].begin(), adj[i].end());

        // aseguramos diagonal al final
        auto it = std::find(cols.begin(), cols.end(), (uint32_t)i);
        if (it == cols.end()) {
            cols.push_back((uint32_t)i);
        } else {
            uint32_t diag = *it;
            cols.erase(it);
            cols.push_back(diag);
        }

        P.row_start[i] = (uint32_t)nnz;
        nnz += cols.size();
    }
    P.row_start[ndof] = (uint32_t)nnz;

    P.col.resize(nnz);
    size_t offset = 0;
    for (size_t i = 0; i < ndof; ++i) {
        auto &cols = rows[i];
        for (size_t k = 0; k < cols.size(); ++k) {
            P.col[offset + k] = cols[k];
        }
        offset += cols.size();
    }

    P.nnz = nnz;
}




void build_P2_mass_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &M)
{
    size_t vtx_count = m.vertex_count();
    size_t tri_count = m.triangle_count();

    // reconstruimos edge numbering igual que en el patrón
    std::map<std::pair<uint32_t,uint32_t>, uint32_t> edge_id;
    uint32_t edge_count = 0;
    build_edge_numbering(m, edge_id, edge_count);

    size_t ndof = vtx_count + edge_count;

    M.symmetric = true;
    M.rows = M.cols = ndof;
    M.nnz = P.nnz;
    M.row_start = P.row_start.data;
    M.col = P.col.data;
    M.data.resize(M.nnz);
    for (size_t i = 0; i < M.nnz; ++i)
        M.data[i] = 0.0;

    const TArray<uint32_t> &idx = m.indices;

    for (size_t t = 0; t < tri_count; ++t) {
        uint32_t a = idx[3*t+0];
        uint32_t b = idx[3*t+1];
        uint32_t c = idx[3*t+2];

        Vec3f A = m.positions[a];
        Vec3f B = m.positions[b];
        Vec3f C = m.positions[c];

        Vec3d AB = { double(B[0]-A[0]), double(B[1]-A[1]), double(B[2]-A[2]) };
        Vec3d AC = { double(C[0]-A[0]), double(C[1]-A[1]), double(C[2]-A[2]) };

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

        for (int i = 0; i < 6; ++i)
            for (int j = 0; j <= i; ++j)
                M(dof[i], dof[j]) += Mloc[i*6 + j];
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
    S.nnz = P.col.size;
    S.row_start = P.row_start.data;
    S.col = P.col.data;
    S.data.resize(S.nnz);
    for (size_t i = 0; i < S.nnz; ++i) {
        S.data[i] = 0.0;
    }

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
