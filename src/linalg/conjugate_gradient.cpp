#include <assert.h>
#include "conjugate_gradient.h"
#include "tiny_blas.h"

double cg_iterate_once(const Matrix &A, double *__restrict x,
		       double *__restrict r, double *__restrict p,
		       double *__restrict Ap, double r2){
	size_t N = A.rows;
	A.mvp(p,Ap); // Ap = A*p
	double alpha = r2/blas_dot(p,Ap,N); //alpha = r2/p2_A
	blas_axpy(alpha,p,x,N); // x = x + alpha*p
	blas_axpy(-alpha,Ap,r,N); // r = r -alpha*A*p
	double new_r2 = blas_dot(r,r,N); //r2_{n+1}
	double beta = new_r2/r2; // beta = r2_{n+1}/r2_n
	blas_axpby(1,r,beta,p,N); // p = r + beta p
	return new_r2;
}


size_t conjugate_gradient_solve(const Matrix &A, const double *__restrict b,
				double *__restrict x, double *__restrict r,
				double *__restrict p, double *__restrict Ap,
				double *__restrict rel_error, double tol,
				int max_iter, bool inited = false){
	assert(A.rows == A.cols);
	size_t N = A.rows;
	
	if (inited){
		// r0 = b - Ax0
		A.mvp(x,r);
		blas_axpby(1,b,-1,r,N);
		blas_copy(r,p,N);
	}
	double b2 = blas_dot(b,b,N);
	double r2 = blas_dot(r,r,N);
	*rel_error = r2/b2;
	size_t iter = 0;
	while(iter<max_iter && *rel_error>tol){
		r2 = cg_iterate_once(A,x,r,p,Ap,r2);
		*rel_error = r2/b2;
		iter++;
	}
	return iter;
}
