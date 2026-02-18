#include "sparse_matrix.h"
#include "array.h"


/* Prerforms Cholesky factorization of a sparse matrix A, provided the pre-computed
sparsity pattern of the factor
Input:
A - SPD matrix
P - Cholesky sparsity pattern associated to A
L - Sparse matrix where result will be stored*/
void cholesky_factorization(CSRMatrix &A, CSRPattern&P,  CSRMatrix &L);

/* Solve the linear system L^tLx = b where L is the Cholesky factor of an SPD matrix A.
Input:
L - lower triangular Cholesky factor
b - right hand side
x - solution vector*/
void solve_cholesky(CSRMatrix &L,  const double *__restrict b, double *__restrict x);
