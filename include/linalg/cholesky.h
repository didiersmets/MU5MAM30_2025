#pragma once

#include "sparse_matrix.h"
#include "array.h"

/******************************************************************************
 * Sparse Cholesky factorization  A = L * L^T
 *
 * L is stored as:
 *   - A CSRPattern  L_pat  that OWNS the row_start[] and col[] arrays.
 *   - A CSRMatrix   L      whose row_start/col RAW POINTERS point into
 *     L_pat's TArrays, and whose data[] holds the numerical values.
 *
 * Only the LOWER triangle is stored (including diagonal).
 * Within each row, column indices are sorted ascending (diagonal last).
 *
 * Usage:
 *   CholeskySolver chol;
 *   chol.factorize(A);     // once before the simulation loop
 *   chol.solve(b, x);      // at every time step, replaces CG
 ******************************************************************************/
struct CholeskySolver {

    size_t N = 0;

    /* Elimination tree (Algorithm 4.2).
     * parent[i] = parent of node i; parent[root] = N (sentinel).         */
    TArray<uint32_t> parent;

    /* Sparsity pattern of L — owns the memory.                            */
    CSRPattern L_pat;

    /* Numerical factor L.
     * L.row_start and L.col are raw pointers into L_pat.row_start.data
     * and L_pat.col.data respectively.  L.data holds the doubles.        */
    CSRMatrix L;

    /* Factorize A = L L^T (symbolic then numeric).
     * A must be SPD and stored lower-triangular with symmetric = true.    */
    void factorize(const CSRMatrix &A);

    /* Solve A x = b via forward + backward substitution with L.           */
    void solve(const double *b, double *x) const;

private:
    void build_etree  (const CSRMatrix &A);  /* Algorithm 4.2 */
    void build_pattern(const CSRMatrix &A);  /* Algorithm 4.3 */
    void numeric      (const CSRMatrix &A);  /* Algorithm 5.6 */

    void forward_solve (const double *b, double *y) const;
    void backward_solve(const double *y, double *x) const;
};