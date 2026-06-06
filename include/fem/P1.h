#pragma once

#include "fem_matrix.h"
#include "sparse_matrix.h"

/* CSR variants */
void build_P1_CSRPattern(const Mesh &m, CSRPattern &P);
void build_P1_mass_matrix(const Mesh &m, const CSRPattern &P, CSRMatrix &M);
void build_P1_stiffness_matrix(const Mesh &m, const CSRPattern &P,
			       CSRMatrix &S);

/* Nonlinear minimal graph helpers (CSR variants) */
void build_P1_stiffness_matrix_NS(const Mesh &m, const CSRPattern &P,
				  CSRMatrix &S, const double *den,
				  const double *u,
				  double area);
void build_P1_rhs_NS(const Mesh &m, const double *den, const double *u,
		     TArray<double> &rhs);

/* FEM matrix variants */
void build_P1_mass_matrix(const Mesh &m, FEMatrix &M);
void build_P1_stiffness_matrix(const Mesh &m, FEMatrix &S);
