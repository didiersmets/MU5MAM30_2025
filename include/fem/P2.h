#pragma once

#include "fem_matrix.h"
#include "sparse_matrix.h"

#include "adjacency.h"
#include "adjacency_edge.h"

/* CSR variants */
void build_P2_CSRPattern_vertex_rows(const Mesh &m, const VTAdjacency &vt_adj, size_t &nnz, CSRPattern &P);
void build_P2_CSRPattern_edge_rows(const Mesh &m, const ETAdjacency &et_adj, size_t &nnz, CSRPattern &P);
void build_P2_CSRPattern(const Mesh &m, CSRPattern &P);

void build_P2_mass_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &M);
void build_P2_stiffness_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &S);
