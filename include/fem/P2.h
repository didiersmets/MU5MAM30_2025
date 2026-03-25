#pragma once

#include <vector>

#include "fem_matrix.h"
#include "sparse_matrix.h"
#include "edge.h"

/* CSR variants */
void build_P2_CSRPattern(const Mesh &m, CSRPattern &P, std::vector<Edge> &out_edges);
void build_P2_mass_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &M, std::vector<Edge> &out_edges);
void build_P2_stiffness_matrix(const Mesh &m, const CSRPattern &P,
			       CSRMatrix &S, std::vector<Edge> &out_edges);