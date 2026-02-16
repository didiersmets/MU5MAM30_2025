#include <stdio.h>
#include <string.h>
#include <iostream>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "P2.h"
#include "fem_matrix.h"
#include "mass_P2.h"
#include "mesh.h"
#include "sparse_matrix.h"
#include "stiffness_P2.h"
#include "hash_table.h"
#include "hash.h"

#include "hashers.h"

#include <unordered_set>

using namespace std;

/* Fill the vertex rows of the sparsity pattern */
void build_P2_CSRPattern_vertex_rows(const Mesh &m, const VTAdjacency &vt_adj, size_t &nnz, CSRPattern &P)
{
    size_t nv = m.vertex_count();

    VertexPairHasher hasher_vv{};
    VertexPair current_vertex_pair;

    uint32_t *current_key;
    uint32_t dummy = 0;

    for (uint32_t k = 0; k < nv; k++)
    {
        HashTable<VertexPair, uint32_t, VertexPairHasher> seen_vv(2 * vt_adj.degree[k] + 4, hasher_vv);

        P.row_start[k] = nnz;

        /* The first vertex is the current one */
        current_vertex_pair.first = k;

        /* Add the diagonal */
        P.col[nnz++] = k;
        seen_vv.set_at({k, k}, dummy);

        for (size_t j = vt_adj.offset[k]; j < vt_adj.offset[k] + vt_adj.degree[k]; j++)
        {
            uint32_t next = vt_adj.vtri[j].next;
            uint32_t prev = vt_adj.vtri[j].prev;

            uint32_t second_vertex[2] = {next, prev};

            for (size_t l = 0; l < 2; l++)
            {
                current_vertex_pair.second = second_vertex[l];
                current_key = seen_vv.get(current_vertex_pair);

                /* Keep only the new interactions, belonging to the lower triangular part of the pattern */
                if (!current_key && current_vertex_pair.second <= k)
                {
                    P.col[nnz++] = current_vertex_pair.second;
                    seen_vv.set_at(current_vertex_pair, dummy);
                }
            }
        }
    }
    P.row_start[nv] = nnz;
}

/* Fill the edge rows of the sparsity pattern */
void build_P2_CSRPattern_edge_rows(const Mesh &m, const ETAdjacency &et_adj, size_t &nnz, CSRPattern &P)
{
    size_t nv = m.vertex_count();
    size_t ne = m.edge_count();

    VertexEdgePairHasher hasher_ve{};
    VertexEdgePair current_vertex_edge_pair;

    EdgePairHasher hasher_ee{};
    EdgePair current_edge_pair;

    uint32_t *current_key;
    uint32_t dummy = 0;

    for (size_t k = 0; k < ne; k++)
    {
        HashTable<VertexEdgePair, uint32_t, VertexEdgePairHasher> seen_ve(2 * et_adj.degree[k] + 4, hasher_ve);
        HashTable<EdgePair, uint32_t, EdgePairHasher> seen_ee(2 * et_adj.degree[k] + 4, hasher_ee);

        P.row_start[k + nv] = nnz;

        /* The edge is the current one */
        current_vertex_edge_pair.edge = m.edges[k];

        /* Add the vertex-edge contributions */
        for (size_t j = et_adj.offset[k]; j < et_adj.offset[k] + et_adj.degree[k]; j++)
        {
            uint32_t tri = et_adj.etri[j].tri;

            uint32_t a = m.indices[3 * tri];
            uint32_t b = m.indices[3 * tri + 1];
            uint32_t c = m.indices[3 * tri + 2];

            uint32_t current_vertex[3] = {a, b, c};

            for (size_t l = 0; l < 3; l++)
            {
                current_vertex_edge_pair.vertex_idx = current_vertex[l];
                current_key = seen_ve.get(current_vertex_edge_pair);

                if (!current_key)
                {
                    P.col[nnz++] = current_vertex[l];
                    seen_ve.set_at(current_vertex_edge_pair, dummy);
                }
            }
        }

        /* The first edge is the current one */
        current_edge_pair.first = m.edges[k];
        current_edge_pair.second = m.edges[k];

        /* Add the diagonal */
        P.col[nnz++] = nv + *m.edge_idx.get(m.edges[k]);
        seen_ee.set_at({current_edge_pair.first, current_edge_pair.second}, dummy);

        /* Add the edge-edge contributions */
        for (size_t j = et_adj.offset[k]; j < et_adj.offset[k] + et_adj.degree[k]; j++)
        {
            /* The second edge is either the current edge, or the two other edges */
            Edge second_edge[2] = {et_adj.etri[j].next, et_adj.etri[j].prev};

            for (size_t l = 0; l < 2; l++)
            {
                current_edge_pair.second = second_edge[l];
                current_key = seen_ee.get(current_edge_pair);

                uint32_t col = nv + *m.edge_idx.get(current_edge_pair.second);
                uint32_t row = nv + k;

                /* Keep only the new interactions, belonging to the lower triangular part of the pattern */
                if (!current_key && col <= row)
                {
                    P.col[nnz++] = col;
                    seen_ee.set_at(current_edge_pair, dummy);
                }
            }
        }
    }
    P.row_start[nv + ne] = nnz;
}

/* CSRMatrix variants */
void build_P2_CSRPattern(const Mesh &m, CSRPattern &P)
{
    /* Adjacency mappings */
    VTAdjacency vt_adj(m);
    ETAdjacency et_adj(m);

    size_t nv = m.vertex_count();
    size_t ne = m.edge_count();

    P.symmetric = true;
    P.rows = nv + ne;
    P.cols = nv + ne;

    size_t nnz = 0;
    size_t nnz_non_symmetric = 2 * (nv + ne);
    for (size_t k = 0; k < nv; k++)
        nnz_non_symmetric += vt_adj.degree[k];

    for (size_t k = 0; k < ne; k++)
        nnz_non_symmetric += et_adj.degree[k];

    P.row_start.resize(nv + ne + 1);
    P.col.resize(2 * nnz_non_symmetric);

    /* Build the vertex rows (containing vertex-vertex interactions) */
    build_P2_CSRPattern_vertex_rows(m, vt_adj, nnz, P);

    /* Build the bloc edge rows (containing edge-vertex and edge-edge interactions) */
    build_P2_CSRPattern_edge_rows(m, et_adj, nnz, P);

    P.col.resize(nnz);
    P.nnz = nnz;

    /* Reorder P.col */
    for (size_t i = 0; i < nv + ne; i++)
    {
        size_t start = P.row_start[i];
        size_t end = P.row_start[i + 1];

        for (size_t j = start; j < end; j++)
        {
            for (size_t k = j + 1; k < end; k++)
            {
                if (P.col[k] < P.col[j])
                {
                    uint32_t temp = P.col[j];
                    P.col[j] = P.col[k];
                    P.col[k] = temp;
                }
            }
        }
    }
}

void build_P2_mass_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &M)
{
    size_t nv = m.vertex_count();
    assert(P.row_start.size == nv + m.edge_count() + 1);

    M.symmetric = P.symmetric;
    M.rows = M.cols = P.rows;
    M.nnz = P.col.size;
    M.row_start = P.row_start.data;
    M.col = P.col.data;
    M.data.resize(M.nnz);
    for (size_t i = 0; i < M.nnz; ++i)
    {
        M.data[i] = 0.0;
    }

    double mass_matrix[36];

    /* Build the global mass matrix */
    size_t nt = m.triangle_count();
    for (size_t tri = 0; tri < nt; tri++)
    {
        uint32_t a = m.indices[3 * tri];
        uint32_t b = m.indices[3 * tri + 1];
        uint32_t c = m.indices[3 * tri + 2];

        Vec3 A = m.positions[a];
        Vec3 B = m.positions[b];
        Vec3 C = m.positions[c];

        Edge e_ab(a, b);
        Edge e_bc(b, c);
        Edge e_ca(c, a);

        uint32_t ab = *m.edge_idx.get(e_ab);
        uint32_t bc = *m.edge_idx.get(e_bc);
        uint32_t ca = *m.edge_idx.get(e_ca);

        Vec3d AB = {(double)B[0] - (double)A[0],
                    (double)B[1] - (double)A[1],
                    (double)B[2] - (double)A[2]};
        Vec3d AC = {(double)C[0] - (double)A[0],
                    (double)C[1] - (double)A[1],
                    (double)C[2] - (double)A[2]};

        mass_P2(AB, AC, mass_matrix);

        /* Bloc vertex-vertex */
        M(a, a) += mass_matrix[0];
        M(b, b) += mass_matrix[7];
        M(c, c) += mass_matrix[14];

        if (b < a)
            M(a, b) += mass_matrix[1];
        else
            M(b, a) += mass_matrix[1];

        if (c < a)
            M(a, c) += mass_matrix[2];
        else
            M(c, a) += mass_matrix[2];

        if (c < b)
            M(b, c) += mass_matrix[8];
        else
            M(c, b) += mass_matrix[8];

        /* Bloc vertex-edge */
        M(nv + ab, a) += mass_matrix[3];
        M(nv + bc, a) += mass_matrix[4];
        M(nv + ca, a) += mass_matrix[5];

        M(nv + ab, b) += mass_matrix[9];
        M(nv + bc, b) += mass_matrix[10];
        M(nv + ca, b) += mass_matrix[11];

        M(nv + ab, c) += mass_matrix[15];
        M(nv + bc, c) += mass_matrix[16];
        M(nv + ca, c) += mass_matrix[17];

        /* Bloc edge-edge */
        M(nv + ab, nv + ab) += mass_matrix[21];
        M(nv + bc, nv + bc) += mass_matrix[28];
        M(nv + ca, nv + ca) += mass_matrix[35];

        if (bc < ab)
            M(nv + ab, nv + bc) += mass_matrix[22];
        else
            M(nv + bc, nv + ab) += mass_matrix[22];

        if (ca < ab)
            M(nv + ab, nv + ca) += mass_matrix[23];
        else
            M(nv + ca, nv + ab) += mass_matrix[23];

        if (ca < bc)
            M(nv + bc, nv + ca) += mass_matrix[29];
        else
            M(nv + ca, nv + bc) += mass_matrix[29];
    }
}

void build_P2_stiffness_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &S)
{
    size_t nv = m.vertex_count();
    assert(P.row_start.size == nv + m.edge_count() + 1);

    S.symmetric = P.symmetric;
    S.rows = S.cols = P.rows;
    S.nnz = P.col.size;
    S.row_start = P.row_start.data;
    S.col = P.col.data;
    S.data.resize(S.nnz);
    for (size_t i = 0; i < S.nnz; ++i)
    {
        S.data[i] = 0.0;
    }

    double stiffness_matrix[36];

    /* Build the global stiffness matrix */
    size_t nt = m.triangle_count();
    for (size_t tri = 0; tri < nt; tri++)
    {
        uint32_t a = m.indices[3 * tri];
        uint32_t b = m.indices[3 * tri + 1];
        uint32_t c = m.indices[3 * tri + 2];

        Vec3 A = m.positions[a];
        Vec3 B = m.positions[b];
        Vec3 C = m.positions[c];

        Edge e_ab(a, b);
        Edge e_bc(b, c);
        Edge e_ca(c, a);

        uint32_t ab = *m.edge_idx.get(e_ab);
        uint32_t bc = *m.edge_idx.get(e_bc);
        uint32_t ca = *m.edge_idx.get(e_ca);

        Vec3d AB = {(double)B[0] - (double)A[0],
                    (double)B[1] - (double)A[1],
                    (double)B[2] - (double)A[2]};
        Vec3d AC = {(double)C[0] - (double)A[0],
                    (double)C[1] - (double)A[1],
                    (double)C[2] - (double)A[2]};

        stiffness_P2(AB, AC, stiffness_matrix);

        /* Bloc vertex-vertex */
        S(a, a) += stiffness_matrix[0];
        S(b, b) += stiffness_matrix[7];
        S(c, c) += stiffness_matrix[14];

        if (b < a)
            S(a, b) += stiffness_matrix[1];
        else
            S(b, a) += stiffness_matrix[1];

        if (c < a)
            S(a, c) += stiffness_matrix[2];
        else
            S(c, a) += stiffness_matrix[2];

        if (c < b)
            S(b, c) += stiffness_matrix[8];
        else
            S(c, b) += stiffness_matrix[8];

        /* Bloc vertex-edge */
        S(nv + ab, a) += stiffness_matrix[3];
        S(nv + bc, a) += stiffness_matrix[4];
        S(nv + ca, a) += stiffness_matrix[5];

        S(nv + ab, b) += stiffness_matrix[9];
        S(nv + bc, b) += stiffness_matrix[10];
        S(nv + ca, b) += stiffness_matrix[11];

        S(nv + ab, c) += stiffness_matrix[15];
        S(nv + bc, c) += stiffness_matrix[16];
        S(nv + ca, c) += stiffness_matrix[17];

        /* Bloc edge-edge */
        S(nv + ab, nv + ab) += stiffness_matrix[21];
        S(nv + bc, nv + bc) += stiffness_matrix[28];
        S(nv + ca, nv + ca) += stiffness_matrix[35];

        if (bc < ab)
            S(nv + ab, nv + bc) += stiffness_matrix[22];
        else
            S(nv + bc, nv + ab) += stiffness_matrix[22];

        if (ca < ab)
            S(nv + ab, nv + ca) += stiffness_matrix[23];
        else
            S(nv + ca, nv + ab) += stiffness_matrix[23];

        if (ca < bc)
            S(nv + bc, nv + ca) += stiffness_matrix[29];
        else
            S(nv + ca, nv + bc) += stiffness_matrix[29];
    }
}
