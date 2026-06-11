#pragma once
#include <stdint.h>
#include "matrix.h"

/* =========================================================================
 * Unpreconditioned CG
 * ========================================================================= */

double cg_iterate_once(const Matrix &A, double *__restrict x,
                       double *__restrict r, double *__restrict p,
                       double *__restrict Ap, double r2);

size_t conjugate_gradient_solve(const Matrix &A, const double *__restrict b,
                                double *__restrict x, double *__restrict r,
                                double *__restrict p, double *__restrict Ap,
                                double *rel_error, double tol,
                                int max_iter, bool inited = false);

/* =========================================================================
 * Preconditioners
 *
 * Templated on the concrete matrix type (e.g. CSRMatrix) so that element
 * access (row_start, col, data arrays) is available at build time.
 * The abstract Matrix base only provides mvp() and cannot be used here.
 *
 * Each preconditioner exposes:
 *   build(const Mat &A)                    -- compute factors once
 *   apply(const double *r, double *z)      -- solve M z = r
 *   free_data()                            -- release owned memory
 * ========================================================================= */

template <typename Mat>
struct JacobiPrecond {
        double  *diag = nullptr;
        size_t   N    = 0;

        void build(const Mat &A);
        void apply(const double *__restrict r, double *__restrict z) const;
        void free_data();
};

template <typename Mat>
struct SSORPrecond {
        const Mat *A    = nullptr;
        double    *diag = nullptr;  /* cached diagonal for O(1) access */
        size_t     N    = 0;

        void build(const Mat &A);
        void apply(const double *__restrict r, double *__restrict z) const;
        void free_data();
};

template <typename Mat>
struct IncompleteCholeskyPrecond {
        uint32_t *row_start = nullptr;  /* same sparsity as upper(A) */
        uint32_t *col       = nullptr;
        double   *data      = nullptr;
        size_t    N         = 0;
        size_t    nnz       = 0;

        void build(const Mat &A);
        void apply(const double *__restrict r, double *__restrict z) const;
        void free_data();
};

/* =========================================================================
 * Preconditioned CG
 *
 * Same signature as conjugate_gradient_solve plus:
 *   z  -- extra scratch vector of length N  (holds M^{-1} r each step)
 *   M  -- a built preconditioner
 * ========================================================================= */
template <typename Precond>
size_t pcg_solve(const Matrix &A, const double *__restrict b,
                 double *__restrict x,  double *__restrict r,
                 double *__restrict z,  double *__restrict p,
                 double *__restrict Ap,
                 double *rel_error, double tol, int max_iter,
                 bool inited, const Precond &M);

/* Explicit instantiation declarations for CSRMatrix (defined in .cpp) */
struct CSRMatrix;

extern template struct JacobiPrecond<CSRMatrix>;
extern template struct SSORPrecond<CSRMatrix>;
extern template struct IncompleteCholeskyPrecond<CSRMatrix>;

extern template size_t pcg_solve<JacobiPrecond<CSRMatrix>>(
        const Matrix &, const double *, double *, double *, double *,
        double *, double *, double *, double, int, bool,
        const JacobiPrecond<CSRMatrix> &);

extern template size_t pcg_solve<SSORPrecond<CSRMatrix>>(
        const Matrix &, const double *, double *, double *, double *,
        double *, double *, double *, double, int, bool,
        const SSORPrecond<CSRMatrix> &);

extern template size_t pcg_solve<IncompleteCholeskyPrecond<CSRMatrix>>(
        const Matrix &, const double *, double *, double *, double *,
        double *, double *, double *, double, int, bool,
        const IncompleteCholeskyPrecond<CSRMatrix> &);