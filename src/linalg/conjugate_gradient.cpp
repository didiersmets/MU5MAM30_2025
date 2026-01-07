#include "linalg/conjugate_gradient.h"

#include "linalg/tiny_blas.h"
#include "matrix/matrix.h"
#include "matrix/sparse_matrix.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

double cg_iterate_once(const Matrix& A,
                       double* __restrict x,
                       double* __restrict r,
                       double* __restrict p,
                       double* __restrict Ap,
                       double r2)
{
  /* We use the following iteration process :
   * \alpha_n = ||r_n||^2 / ||p_n||_A^2
   * x_{n+1} =  x_n + \alpha_n * p_n
   * r_{n+1} = r_n - \alpha_n * A * p_n
   * \beta_{n+1} = ||r_{n+1}||^2 / ||r_n||^2
   * p_{n+1} = r_{n+1} + \beta_{n+1} * p_n
   */

  size_t N = A.rows;

  /* Computation of A * p_n */
  A.mvp(p, Ap);

  /* Computation of \alpha_n */
  double p_A2  = blas_dot(p, Ap, N);
  double alpha = r2 / p_A2;

  /* Computation of x_{n+1} */
  blas_axpby(alpha, p, 1., x, N);

  /* Computation of r_{n+1} */
  blas_axpby(-alpha, Ap, 1., r, N);

  /* Computation of \beta_{n+1} */
  double new_r2 = blas_dot(r, r, N);
  double beta   = new_r2 / r2;

  /* Computation of p_{n+1} */
  blas_axpby(1., r, beta, p, N);

  return new_r2;
}

size_t conjugate_gradient_solve(const Matrix& A,
                                const double* __restrict b,
                                double* __restrict x,
                                double* __restrict r,
                                double* __restrict p,
                                double* __restrict Ap,
                                double* rel_error,
                                double  tol,
                                int     max_iter,
                                bool    inited)
{
  size_t N = A.rows;
  assert(A.rows == A.cols);

  double b2 = blas_dot(b, b, N);

  if (!inited)
  {
    /* r_0 = b - Ax_0 */
    A.mvp(x, r);
    blas_axpby(1, b, -1, r, N);
    /* p_0 = r_0 */
    blas_copy(r, p, N);
  }

  double r2  = blas_dot(r, r, N);
  *rel_error = sqrt(r2 / b2);

  int iter = 0;
  while ((iter < max_iter) && (*rel_error > tol))
  {
    r2         = cg_iterate_once(A, x, r, p, Ap, r2);
    *rel_error = sqrt(r2 / b2);
    iter++;
  }
  return iter;
}
