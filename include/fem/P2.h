#pragma once

#include "mesh.h"
#include "sparse_matrix.h"

#include <unordered_map>

/* CSR variants */
void build_P2_CSRPattern(const Mesh &m,
			 CSRPattern &P,
			 std::unordered_map<uint64_t, uint32_t> &edge2dof);
void build_P2_mass_matrix(const Mesh &m,
			  const CSRPattern &P,
			  CSRMatrix &M,
			  std::unordered_map<uint64_t, uint32_t> &edge2dof);
void build_P2_stiffness_matrix(const Mesh &m,
			       const CSRPattern &P,
			       CSRMatrix &S,
			       std::unordered_map<uint64_t, uint32_t> &edge2dof);
