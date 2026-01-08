#include "matrix/matrix.h"

/* The conjugate gradient algorithm to solve A * x = b for a SPD square
 * matrix A of size NxN rely on a number of vectors and one matrix
 * vector computation.
 * In the implementation below we let the user responsible
 * for allocation of memory for the corresponding vectors, including
 * the internal ones : x is the solution vector, b the rhd, r is the
 * residual vector (r = b - A * x), p is the direction of search
 * vector, and Ap to store the product A*p.
 */

/* Make one iteration (i.e. update x, r, p and Ap, and return the updated
 * value for the (squared, hence the name r2) norm or the residual.
 */
double cg_iterate_once(const Matrix& A,
                       double* __restrict x,
                       double* __restrict r,
                       double* __restrict p,
                       double* __restrict Ap,
                       double r2);

/* The master function, that call cg_interate_once as long as "needed",
 * the latter being guided by a user defined max iteration count, and
 * tolerance to error. Inited = false correspond to the whole classical
 * CG algorithm, inited = true alows to start from a given iterate (the value
 * of r and p are not reset from b, but considered known), this allows e.g.
 * to observe the algorithm by batches of a few iterates.
 */
size_t conjugate_gradient_solve(const Matrix& A,
                                const double* __restrict b,
                                double* __restrict x,
                                double* __restrict r,
                                double* __restrict p,
                                double* __restrict Ap,
                                double* __restrict rel_error,
                                double tol,
                                int    max_iter,
                                bool   inited = false);
