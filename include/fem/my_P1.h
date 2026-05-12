#pragma once

#include "my_sparse_matrix.h"
#include "my_mesh.h"
#include "array.h"

void build_P1_CSRPattern(const MyMesh &m, MyCSRPattern &P);
void build_P1_stiffness_matrix(const MyMesh &m, const MyCSRPattern &P,
			       MyCSRMatrix &S, bool modified = false, const double *den = nullptr);
void build_P1_stiffness_matrix_NS(const MyMesh &m, const MyCSRPattern &P,
				   MyCSRMatrix &S, const double *den, const double *u);
void build_P1_rhs_NS(const MyMesh &m, const double *den, const double *u,
				   TArray<double> &rhs);