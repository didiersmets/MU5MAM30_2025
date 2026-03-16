#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <vector>

#include "projected_conjugate_gradient.h"
#include "matrix.h"
#include "tiny_blas.h"

// this methode is for boundary conditions 
static inline void project_boundary(double *vec, size_t N, bool has_boundary, const std::vector<bool> *is_boundary) 
{
    if (has_boundary && is_boundary != nullptr) {
        for (size_t i = 0; i < N; ++i) {
            if ((*is_boundary)[i]) {
                vec[i] = 0.0;
            }
        }
    }
}

double projected_cg_iterate_once(const Matrix &A, double *__restrict x,
                                 double *__restrict r, double *__restrict p,
                                 double *__restrict Ap, double r2,
                                 bool has_boundary, const std::vector<bool> *is_boundary)
{
    size_t N = A.rows;

    A.mvp(p, Ap);
    double pAp = blas_dot(p, Ap, N);
    double alpha = r2 / pAp;

    blas_axpy(alpha, p, x, N);
    blas_axpy(-alpha, Ap, r, N);

    // projection
    project_boundary(r, N, has_boundary, is_boundary);

    double r2_new = blas_dot(r, r, N);
    double beta = r2_new / r2;
    blas_axpby(1.0, r, beta, p, N);

    project_boundary(p, N, has_boundary, is_boundary);

    return r2_new;
}

size_t projected_conjugate_gradient_solve(const Matrix &A, const double *__restrict b,
                                          double *__restrict x, double *__restrict r,
                                          double *__restrict p, double *__restrict Ap,
                                          double *rel_error, double tol, int max_iter,
                                          bool inited,
                                          bool has_boundary, const std::vector<bool> *is_boundary)
{
    size_t N = A.rows;
    assert(A.rows == A.cols);

    double b2 = blas_dot(b, b, N);

    if (!inited) {
        A.mvp(x, r);
        blas_axpby(1, b, -1, r, N);

        project_boundary(r, N, has_boundary, is_boundary);

        blas_copy(r, p, N);
    }

    double r2 = blas_dot(r, r, N);
    
    if (b2 == 0.0) {
        *rel_error = 0.0;
        return 0;
    }
    
    *rel_error = sqrt(r2 / b2);

    size_t iter = 0;
    while ((iter < (size_t)max_iter) && (*rel_error > tol)) {
        r2 = projected_cg_iterate_once(A, x, r, p, Ap, r2, has_boundary, is_boundary);
        *rel_error = sqrt(r2 / b2);
        iter++;
    }
    
    return iter;
}