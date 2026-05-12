#pragma once

#include "matrix.h"
#include "my_sparse_matrix.h"

double cg_iterate_once(const MyCSRMatrix &A, double *__restrict x,
		       double *__restrict r, double *__restrict p,
		       double *__restrict Ap, double r2);

size_t conjugate_gradient_solve(const MyCSRMatrix &A, const double *__restrict b,
				double *__restrict x, double *__restrict r,
				double *__restrict p, double *__restrict Ap,
				double *__restrict rel_error, double tol,
				int max_iter, bool inited = false);
