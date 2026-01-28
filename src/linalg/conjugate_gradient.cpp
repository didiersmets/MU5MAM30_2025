#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "conjugate_gradient.h"
#include "matrix.h"
#include "tiny_blas.h"

double cg_iterate_once(const Matrix &A,
					   double *__restrict x,
					   double *__restrict r,
					   double *__restrict p,
					   double *__restrict Ap,
					   double r2)
{
	const size_t N = A.rows;
	A.mvp(p, Ap);

	//alpha k-1
	double pAp = blas_dot(p, Ap, N);//<p_{k-1}, Ap_{k-1}>
	double alpha = r2 / pAp;

	// mise à jour de x => x<- x + alpha*p (on a une forme axpy )
	blas_axpy(alpha, p, x, N);

	// mise à jour de r => r<- r - alpha Ap
	blas_axpy(-alpha, Ap, r, N);

	//on calcul le nouveau r², on a toujours le ||r_{k-1}||**2 en stock dans r2
	double r2_new = blas_dot(r, r, N);

	// on calcul beta qui est  ||r_k||**2 / ||r_{k-1}||**2
	double beta = r2_new / r2;

	// On met à jorus pk => p<-r + beta*p
	blas_axpby(1.0, r, beta, p, N);




	return r2_new;
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
