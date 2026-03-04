#pragma once

#include "mesh.h"
#include "sparse_matrix.h"

/* CSR variants */
void build_P2_CSRPattern(Mesh &m,
			 CSRPattern &P);
void build_P2_mass_matrix(const Mesh &m,
			  const CSRPattern &P,
			  CSRMatrix &M);
void build_P2_stiffness_matrix(const Mesh &m,
			       const CSRPattern &P,
			       CSRMatrix &S);
