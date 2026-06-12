#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "conjugate_gradient.h"
#include "matrix.h"
#include "tiny_blas.h"
#include "sparse_matrix.h"

/* =========================================================================
 * Unpreconditioned CG (unchanged)
 * ========================================================================= */

double cg_iterate_once(const Matrix &A, double *__restrict x,
                       double *__restrict r, double *__restrict p,
                       double *__restrict Ap, double r2)
{
        size_t N = A.rows;
        assert(A.rows == A.cols);

        A.mvp(p, Ap);
        double alpha = r2 / blas_dot(p, Ap, N);

        blas_axpy(alpha, p, x, N);
        blas_axpy(-alpha, Ap, r, N);

        double r2_new = blas_dot(r, r, N);
        double beta = r2_new / r2;
        blas_axpby(1, r, beta, p, N);

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
                A.mvp(x, r);
                blas_axpby(1, b, -1, r, N);
                blas_copy(r, p, N);
        }

        double r2 = blas_dot(r, r, N);
        *rel_error = (b2 > 0.0) ? sqrt(r2 / b2) : sqrt(r2);

        int iter = 0;
        while ((iter < max_iter) && (*rel_error > tol)) {
                r2 = cg_iterate_once(A, x, r, p, Ap, r2);
                *rel_error = (b2 > 0.0) ? sqrt(r2 / b2) : sqrt(r2);
                iter++;
        }
        return iter;
}

/* =========================================================================
 * Preconditioners
 *
 * Templated on ConcreteMatrix so that element access (operator(), row_start,
 * col, data) is available.  The abstract Matrix base only provides mvp().
 * ========================================================================= */

/* -------------------------------------------------------------------------
 * Jacobi:  M = diag(A),   z_i = r_i / a_ii
 * ------------------------------------------------------------------------- */
template <typename Mat>
void JacobiPrecond<Mat>::build(const Mat &A)
{
        N    = A.rows;
        diag = (double *)malloc(N * sizeof(double));
        /* Walk the CSR rows to find each diagonal entry */
        for (size_t i = 0; i < N; i++) {
                diag[i] = 0.0;
                for (uint32_t k = A.row_start[i]; k < A.row_start[i+1]; k++) {
                        if (A.col[k] == (uint32_t)i) {
                                diag[i] = A.data[k];
                                break;
                        }
                }
                assert(diag[i] != 0.0);
        }
}

template <typename Mat>
void JacobiPrecond<Mat>::apply(const double *__restrict r,
                               double *__restrict z) const
{
        for (size_t i = 0; i < N; i++)
                z[i] = r[i] / diag[i];
}

template <typename Mat>
void JacobiPrecond<Mat>::free_data()
{
        free(diag);
        diag = nullptr;
}

/* -------------------------------------------------------------------------
 * SSOR (omega=1):  M = (D+L) D^{-1} (D+L^T)
 *
 * Forward sweep  (D+L) t = r  uses only j < i  (lower triangle + diagonal).
 * Backward sweep (D+L^T) z = D t  uses only j > i  (upper triangle).
 * Both sweeps follow CSR row order; the symmetric flag means the stored
 * upper triangle also implicitly represents the lower triangle.
 * ------------------------------------------------------------------------- */
template <typename Mat>
void SSORPrecond<Mat>::build(const Mat &mat)
{
        A = &mat;
        N = mat.rows;
        /* Cache diagonal for O(1) access during sweeps */
        diag = (double *)malloc(N * sizeof(double));
        for (size_t i = 0; i < N; i++) {
                diag[i] = 0.0;
                for (uint32_t k = mat.row_start[i]; k < mat.row_start[i+1]; k++) {
                        if (mat.col[k] == (uint32_t)i) {
                                diag[i] = mat.data[k];
                                break;
                        }
                }
        }
}

template <typename Mat>
void SSORPrecond<Mat>::apply(const double *__restrict r,
                             double *__restrict z) const
{
        double *t = (double *)malloc(N * sizeof(double));

        /* Forward solve (D + L) t = r
         * For a symmetric CSR storing the upper triangle only, the lower
         * triangle entry (i,j) with j<i is the same as stored entry (j,i). */
        for (size_t i = 0; i < N; i++) {
                double s = r[i];
                /* subtract lower-triangle contributions a_ij * t_j, j < i */
                for (uint32_t k = A->row_start[i]; k < A->row_start[i+1]; k++) {
                        uint32_t j = A->col[k];
                        if (j < (uint32_t)i)
                                s -= A->data[k] * t[j];
                }
                t[i] = s / diag[i];
        }

        /* Diagonal scaling: t_i *= a_ii */
        for (size_t i = 0; i < N; i++)
                t[i] *= diag[i];

        /* Backward solve (D + L^T) z = t */
        for (size_t i = N; i-- > 0; ) {
                double s = t[i];
                for (uint32_t k = A->row_start[i]; k < A->row_start[i+1]; k++) {
                        uint32_t j = A->col[k];
                        if (j > (uint32_t)i)
                                s -= A->data[k] * z[j];
                }
                z[i] = s / diag[i];
        }

        free(t);
}

template <typename Mat>
void SSORPrecond<Mat>::free_data()
{
        free(diag);
        diag = nullptr;
        A = nullptr;
}

/* -------------------------------------------------------------------------
 * Incomplete Cholesky IC(0):  M = L L^T
 *
 * Storage layout (from build_P1_CSRPattern / sparse_matrix.cpp):
 *   - LOWER triangle is stored, diagonal is the LAST entry in each row.
 *   - Row i contains columns j < i (strict lower triangle) followed by i.
 *   - Columns within each row are sorted ascending (P1.cpp insertion sort).
 *
 * We compute L such that L L^T ≈ A, where L is lower-triangular.
 * This is the standard "left-looking" IC(0) on the lower triangle.
 *
 * Algorithm (column-based, working in-place on a copy of A's lower triangle):
 *   for k = 0..N-1:
 *     L[k][k] = sqrt(A[k][k])                         (pivot)
 *     for i > k with (i,k) stored:
 *       L[i][k] /= L[k][k]                            (scale column k)
 *     for i > k with (i,k) stored:
 *       for j >= k with (i,j) stored and j <= i:
 *         A[i][j] -= L[i][k] * L[j][k]               (rank-1 update, IC(0): skip if not stored)
 *
 * apply: solve M z = r, i.e. L L^T z = r
 *   Forward:  L t = r   (lower-triangular, all entries in CSR row)
 *   Backward: L^T z = t (upper-triangular; L^T entry (i,j) = L[j][i],
 *                         stored as entry (j,i) in row j of our lower-CSR)
 * ------------------------------------------------------------------------- */
template <typename Mat>
void IncompleteCholeskyPrecond<Mat>::build(const Mat &A)
{
        N   = A.rows;
        nnz = A.nnz;

        row_start = (uint32_t *)malloc((N+1) * sizeof(uint32_t));
        col       = (uint32_t *)malloc(nnz   * sizeof(uint32_t));
        data      = (double   *)malloc(nnz   * sizeof(double));

        memcpy(row_start, A.row_start,      (N+1) * sizeof(uint32_t));
        memcpy(col,       A.col,            nnz   * sizeof(uint32_t));
        memcpy(data,      A.data.data,      nnz   * sizeof(double));

        /* Diagonal is the last entry in each row (guaranteed by CSRPattern).
         * Precompute diagonal indices for O(1) access. */
        uint32_t *diag_idx = (uint32_t *)malloc(N * sizeof(uint32_t));
        for (size_t i = 0; i < N; i++)
                diag_idx[i] = row_start[i+1] - 1;   /* diagonal = last entry */

        const double shift_factor = 1e-8;

        for (size_t k = 0; k < N; k++) {
                /* --- pivot: L[k][k] = sqrt(data[diag_idx[k]]) --- */
                double lkk = data[diag_idx[k]];
                if (lkk <= 0.0)
                        lkk = A.data.data[diag_idx[k]] * shift_factor;
                lkk = sqrt(lkk);
                data[diag_idx[k]] = lkk;

                /* --- scale column k: for each i > k with (i,k) stored,
                 *     L[i][k] /= L[k][k]
                 *     Entry (i,k) lives in row i at some position idx where col[idx]==k. --- */
                for (size_t i = k+1; i < N; i++) {
                        /* scan row i for column k (stops at diagonal = last entry) */
                        for (uint32_t idx = row_start[i]; idx < diag_idx[i]; idx++) {
                                if (col[idx] == (uint32_t)k) {
                                        data[idx] /= lkk;
                                        break;
                                }
                        }
                }

                /* --- rank-1 update: for i > k with (i,k) stored,
                 *     for j in [k+1..i] with (i,j) stored:
                 *       A[i][j] -= L[i][k] * L[j][k]
                 *     IC(0) rule: skip (i,j) pairs not in the stored pattern. --- */
                for (size_t i = k+1; i < N; i++) {
                        /* find L[i][k] in row i */
                        double lik = 0.0;
                        for (uint32_t idx = row_start[i]; idx < diag_idx[i]; idx++) {
                                if (col[idx] == (uint32_t)k) { lik = data[idx]; break; }
                        }
                        if (lik == 0.0) continue;

                        /* update (i,j) for all j in row i with k < j <= i */
                        for (uint32_t idx_j = row_start[i]; idx_j <= diag_idx[i]; idx_j++) {
                                uint32_t j = col[idx_j];
                                if (j <= (uint32_t)k) continue;  /* only j > k */

                                /* find L[j][k] in row j */
                                double ljk = 0.0;
                                for (uint32_t s = row_start[j]; s < diag_idx[j]; s++) {
                                        if (col[s] == (uint32_t)k) { ljk = data[s]; break; }
                                }
                                if (j == (uint32_t)i) ljk = lik;   /* diagonal of update */

                                data[idx_j] -= lik * ljk;
                        }
                }
        }

        free(diag_idx);
}

template <typename Mat>
void IncompleteCholeskyPrecond<Mat>::apply(const double *__restrict r,
                                           double *__restrict z) const
{
        double *t = (double *)malloc(N * sizeof(double));

        /* Forward solve L t = r
         * L is lower-triangular; row i contains columns j < i then diagonal i (last).
         * t[i] = (r[i] - sum_{j<i, (i,j) stored} L[i][j] * t[j]) / L[i][i]  */
        for (size_t i = 0; i < N; i++) {
                double s = r[i];
                uint32_t diag = row_start[i+1] - 1;   /* diagonal is last */
                for (uint32_t idx = row_start[i]; idx < diag; idx++)
                        s -= data[idx] * t[col[idx]];
                t[i] = s / data[diag];
        }

        /* Backward solve L^T z = t
         * L^T is upper-triangular: (L^T)[i][j] = L[j][i] for j > i.
         * Process in reverse: z[j] = (t[j] - sum_{i>j} L[i][j]*z[i]) / L[j][j].
         * Entry L[i][j] (i>j) is stored in row i. We use a reverse scatter:
         * process j = N-1 down to 0, divide z[j] by L[j][j], then subtract
         * L[j][col[idx]] * z[j] from all z[col[idx]] (col[idx] < j). */
        memcpy(z, t, N * sizeof(double));
        for (size_t j = N; j-- > 0; ) {
                uint32_t diag = row_start[j+1] - 1;
                z[j] /= data[diag];
                for (uint32_t idx = row_start[j]; idx < diag; idx++)
                        z[col[idx]] -= data[idx] * z[j];
        }

        free(t);
}

template <typename Mat>
void IncompleteCholeskyPrecond<Mat>::free_data()
{
        free(row_start); row_start = nullptr;
        free(col);       col       = nullptr;
        free(data);      data      = nullptr;
}

/* =========================================================================
 * PCG core
 * ========================================================================= */

template <typename Precond>
static double pcg_iterate_once(const Matrix &A,
                                double *__restrict x,
                                double *__restrict r,
                                double *__restrict z,
                                double *__restrict p,
                                double *__restrict Ap,
                                double rz,
                                const Precond &M)
{
        size_t N = A.rows;
        A.mvp(p, Ap);
        double alpha = rz / blas_dot(p, Ap, N);
        blas_axpy(alpha,  p,  x, N);
        blas_axpy(-alpha, Ap, r, N);
        M.apply(r, z);
        double rz_new = blas_dot(r, z, N);
        double beta   = rz_new / rz;
        blas_axpby(1.0, z, beta, p, N);
        return rz_new;
}

template <typename Precond>
size_t pcg_solve(const Matrix &A, const double *__restrict b,
                 double *__restrict x,  double *__restrict r,
                 double *__restrict z,  double *__restrict p,
                 double *__restrict Ap,
                 double *rel_error, double tol, int max_iter,
                 bool inited, const Precond &M)
{
        size_t N = A.rows;
        assert(A.rows == A.cols);

        double b2 = blas_dot(b, b, N);

        if (!inited) {
                A.mvp(x, r);
                blas_axpby(1.0, b, -1.0, r, N);
                M.apply(r, z);
                blas_copy(z, p, N);
        }

        double rz = blas_dot(r, z, N);
        assert(rz > 0.0);

        double r2  = blas_dot(r, r, N);
        *rel_error = (b2 > 0.0) ? sqrt(r2 / b2) : sqrt(r2);

        int iter = 0;
        while ((iter < max_iter) && (*rel_error > tol)) {
                rz = pcg_iterate_once(A, x, r, z, p, Ap, rz, M);
                r2 = blas_dot(r, r, N);
                *rel_error = (b2 > 0.0) ? sqrt(r2 / b2) : sqrt(r2);
                iter++;
        }
        return iter;
}

/* =========================================================================
 * Explicit instantiations for CSRMatrix
 * ========================================================================= */

template struct JacobiPrecond<CSRMatrix>;
template struct SSORPrecond<CSRMatrix>;
template struct IncompleteCholeskyPrecond<CSRMatrix>;

template size_t pcg_solve<JacobiPrecond<CSRMatrix>>(
        const Matrix &, const double *, double *, double *, double *,
        double *, double *, double *, double, int, bool,
        const JacobiPrecond<CSRMatrix> &);

template size_t pcg_solve<SSORPrecond<CSRMatrix>>(
        const Matrix &, const double *, double *, double *, double *,
        double *, double *, double *, double, int, bool,
        const SSORPrecond<CSRMatrix> &);

template size_t pcg_solve<IncompleteCholeskyPrecond<CSRMatrix>>(
        const Matrix &, const double *, double *, double *, double *,
        double *, double *, double *, double, int, bool,
        const IncompleteCholeskyPrecond<CSRMatrix> &);