#pragma once

#include "matrix/sparse_matrix.h"
#include "mesh/mesh.h"

void build_P1_CSRPattern(const Mesh& m, CSRPattern& P);
void build_P1_mass_matrix(const Mesh& m, const CSRPattern& P, CSRMatrix& M);
void build_P1_stiffness_matrix(const Mesh& m, const CSRPattern& P, CSRMatrix& S);
