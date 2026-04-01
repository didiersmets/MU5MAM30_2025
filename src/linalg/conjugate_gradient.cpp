#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "conjugate_gradient.h"
#include "matrix.h"
#include "tiny_blas.h"

double cg_iterate_once(const Matrix &A, double *__restrict x,
		       double *__restrict r, double *__restrict p,
		       double *__restrict Ap, double r2)
{
	/* Check usage of tiny_blas (BLAS = Basic Linear Algebra Subroutines) 
	 * in include/linalg/tiny_blas.h
	 * BLAS routines are key in efficient linear algebra computations;
	 * we do not really care about best performance here, top notch 
	 * implementations heavilly make use of specific hardware capabilities 
	 * (instruction sets, parallelism, cache efficiency etc). 
	 */

	/* Ap = A*p */
	A.mvp(p, Ap);

	int N = A.rows;
	double alpha = r2 / blas_dot(Ap, p, N);

	/* x_{n+1} = x_n+alpha*p */
	blas_axpy(alpha, p, x, N);

	/* r_{n+1} = r_n - alpha*Ap */
	blas_axpy(-alpha, Ap, r, N);

	double r2_new = blas_dot(r, r, N);
	double beta = r2_new / r2;

	/* p_{n+1} = r_{n+1} + beta*p_n */
	blas_axpby(1, r, beta, p, N);

	return(r2_new);

}

size_t conjugate_gradient_solve(const Matrix &A, const double *__restrict b,
				double *__restrict x, double *__restrict r,
				double *__restrict p, double *__restrict Ap,
				double *rel_error, double tol, int max_iter,
				bool inited)
{
	size_t N = A.rows;
	assert(A.rows == A.cols);

	double b2 = blas_dot(b, b, N);

	if (!inited) {
		/* x_0 = 0 */
    memset(x, 0, N * sizeof(double));
		/* r_0 = b - Ax_0 */
		A.mvp(x, r);
		blas_axpby(1, b, -1, r, N);
		/* p_0 = r_0 */
		blas_copy(r, p, N);
	}

	double r2 = blas_dot(r, r, N);
	*rel_error = sqrt(r2 / b2);

	int iter = 0;
	while ((iter < max_iter) && (*rel_error > tol)) {
		r2 = cg_iterate_once(A, x, r, p, Ap, r2);
		*rel_error = sqrt(r2 / b2);
		iter++;
	}
	return iter;
}
