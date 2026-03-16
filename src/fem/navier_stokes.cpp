#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <cstdio>
#include <cmath>

#include "navier_stokes.h"
#include <limits>

static const uint32_t NIL = std::numeric_limits<uint32_t>::max();

#if !USE_FEM_MATRIX
#include "sparse_cholesky.h"
#endif

#include "P1.h"
#include "tiny_blas.h"



#if !USE_FEM_MATRIX
static void debug_symbolic_phase(const CSRPattern& P)
{
    std::vector<uint32_t> parent;
    CSRPattern PL;

    elimination_tree_from_csr_symmetric(P, parent);

    size_t bad_edges = 0;

    for (uint32_t i = 0; i < P.rows; ++i) {

        size_t start = P.row_start[i];
        size_t stop  = P.row_start[i + 1];

        if (stop == start)
            continue;

        // ignore diagonal (last entry)
        stop--;

        for (size_t k = start; k < stop; ++k) {

            uint32_t j = P.col[k];
            if (j >= i) continue;

            // climb the elimination tree from j
            uint32_t a = j;

            while (a != NIL && a < i)
                a = parent[a];

            if (a != i) {
                printf("ERROR: edge (%u,%u) not represented in elimination tree\n", i, j);
                bad_edges++;
            }
        }
    }

    printf("bad_edges = %zu\n", bad_edges);

    symbolic_cholesky_row_pattern(P, parent, PL);

    size_t nnzA = P.col.size;   // pattern de A (triangle inférieur + diag)
    size_t nnzL = PL.col.size;  // pattern de L (triangle inférieur + diag)
    size_t fill = (nnzL >= nnzA) ? (nnzL - nnzA) : 0;

    printf("\n");
    printf("========== Symbolic phase ==========\n");
    printf("n            = %zu\n", P.rows);
    printf("nnz(A)       = %zu\n", nnzA);
    printf("nnz(L)       = %zu\n", nnzL);
    printf("fill-in      = %zu\n", fill);
    printf("fill ratio   = %.6f\n", nnzA ? (double)nnzL / (double)nnzA : 0.0);

    printf("====================================\n");
    printf("\n");

    for (size_t i = 0; i < PL.rows; ++i) {
        for (size_t k = PL.row_start[i]; k < PL.row_start[i + 1] - 1; ++k) {
            if (PL.col[k] >= i) {
                printf("ERROR: PL has invalid column %u in row %zu\n", PL.col[k], i);
            }
        }
    }
    size_t bad_cols = 0;
    for (size_t i = 0; i < PL.rows; ++i) {
        size_t start = (size_t)PL.row_start[i];
        size_t stop  = (size_t)PL.row_start[i + 1];
        if (stop == start) continue;

        // diag should be last
        if (PL.col[stop - 1] != i) {
            printf("ERROR: diagonal not last in row %zu\n", i);
        }

        // all off-diagonal entries should satisfy j < i
        for (size_t k = start; k + 1 < stop; ++k) {
            if (PL.col[k] >= i) {
                ++bad_cols;
                printf("ERROR: invalid col %u in row %zu\n", PL.col[k], i);
            }
        }

        // sorted row check
        for (size_t k = start + 1; k < stop; ++k) {
            if (PL.col[k - 1] >= PL.col[k]) {
                printf("ERROR: row %zu not strictly increasing at positions %zu/%zu\n",
                    i, k - 1, k);
            }
        }
    }
    printf("bad_cols = %zu\n", bad_cols);

}
#endif

NavierStokesSolver::NavierStokesSolver(const Mesh &m)
	: m(m)
	, N(m.vertex_count())
	, omega(N)
	, Momega(N)
	, psi(N)
#if !USE_FEM_MATRIX
    , ysolve(N) 
#endif
	, r(N)
	, p(N)
	, Ap(N)
{
#if USE_FEM_MATRIX
	build_P1_mass_matrix(m, M);
	build_P1_stiffness_matrix(m, S);
#else
	build_P1_CSRPattern(m, P);
    debug_symbolic_phase(P);
	build_P1_mass_matrix(m, P, M);
	build_P1_stiffness_matrix(m, P, S);
#endif
	vol = M.sum();
	inited = false;
	t = 0;
}

void NavierStokesSolver::set_zero_mean(double *V)
{
	/* Your implementation goes here */
	// Buffer temporaire pour stocker M * V
    std::vector<double> tmp(N);

    M.mvp(V, tmp.data());

    // Approximation de l'intégrale de V
    double integral = blas_sum_in_place(tmp.data(), N);

    // Calcul de la moyenne pondérée
    double mean = integral / vol;

    // Pour avoir intégrale de V = 0
    for (size_t i = 0; i < N; ++i) {
        V[i] -= mean;
    }
}

void NavierStokesSolver::compute_transport(double *T)
{
	memset(T, 0, N * sizeof(double));

	/* Your implementation goes here */
	for (size_t k = 0; k < m.triangle_count(); ++k) {

        // Indices des sommets du triangle k
        const uint32_t ia = m.indices[3 * k + 0];
        const uint32_t ib = m.indices[3 * k + 1];
        const uint32_t ic = m.indices[3 * k + 2];

        assert(ia < N && ib < N && ic < N);

        // Valeurs locales de omega
        const double omega_sum = omega[ia] + omega[ib] + omega[ic];

        // Différences locales de psi
        const double dpsi_bc = psi[ib] - psi[ic];
        const double dpsi_ca = psi[ic] - psi[ia];
        const double dpsi_ab = psi[ia] - psi[ib];

        // Assemblage
        T[ia] += omega_sum * dpsi_bc;
        T[ib] += omega_sum * dpsi_ca;
        T[ic] += omega_sum * dpsi_ab;
    }

    // Normalisation globale
    const double factor = 1.0 / 6.0;
    for (size_t i = 0; i < N; ++i) {
        T[i] *= factor;
    }
}

size_t NavierStokesSolver::compute_stream_function()
{
	size_t iter = 0;

	/* Your implementation goes here */
    double *residual = r.data;
    double *direction = p.data;
    double *Adir = Ap.data;
    double *omega_vec = omega.data;
    double *Momega_vec = Momega.data;
    double *psi_vec = psi.data;

    // b = M * omega
    M.mvp(omega_vec, Momega_vec);

    double rhs_norm2 = blas_dot(Momega_vec, Momega_vec, N);

    // r = b - S * psi
    S.mvp(psi_vec, residual);
    blas_axpby(1.0, Momega_vec, -1.0, residual, N);

    // Direction initiale : p = r
    blas_copy(residual, direction, N);

    double r_norm2 = blas_dot(residual, residual, N);
    double rel_error = sqrt(r_norm2 / rhs_norm2);

    // Boucle CG
    while (rel_error > tol && iter < iter_max) {

        // Ap = S * p
        S.mvp(direction, Adir);

        double denom = blas_dot(direction, Adir, N);
        double alpha = r_norm2 / denom;

        // psi = psi + alpha * p
        blas_axpy(alpha, direction, psi_vec, N);

        // r = r - alpha * Ap
        blas_axpy(-alpha, Adir, residual, N);

        double new_r_norm2 = blas_dot(residual, residual, N);
        rel_error = sqrt(new_r_norm2 / rhs_norm2);

        if (rel_error <= tol)
            break;

        double beta = new_r_norm2 / r_norm2;

        // p = r + beta * p
        blas_axpby(1.0, residual, beta, direction, N);

        r_norm2 = new_r_norm2;
        iter++;
    }


	return iter;
}

#if !USE_FEM_MATRIX
static void debug_check_cholesky_factorization(const CSRMatrix& A, const CSRMatrix& L)
{
    const uint32_t n = (uint32_t)A.rows;

    // Reconstruct dense A and dense LL^T
    std::vector<double> Adense((size_t)n * (size_t)n, 0.0);
    std::vector<double> LLt((size_t)n * (size_t)n, 0.0);

    // Build dense A from lower-triangular symmetric CSR
    for (uint32_t i = 0; i < n; ++i) {
        size_t start = (size_t)A.row_start[i];
        size_t stop  = (size_t)A.row_start[i + 1];
        for (size_t k = start; k < stop; ++k) {
            uint32_t j = A.col[k];
            double aij = A.data[k];

            Adense[(size_t)i * n + j] = aij;
            Adense[(size_t)j * n + i] = aij; 
        }
    }

    // Compute LL^T explicitly in dense form
    for (uint32_t i = 0; i < n; ++i) {
        size_t Li_start = (size_t)L.row_start[i];
        size_t Li_stop  = (size_t)L.row_start[i + 1];

        for (uint32_t j = 0; j < n; ++j) {
            size_t Lj_start = (size_t)L.row_start[j];
            size_t Lj_stop  = (size_t)L.row_start[j + 1];

            double s = 0.0;

            size_t pi = Li_start;
            size_t pj = Lj_start;
            while (pi < Li_stop && pj < Lj_stop) {
                uint32_t ci = L.col[pi];
                uint32_t cj = L.col[pj];

                if (ci == cj) {
                    s += L.data[pi] * L.data[pj];
                    ++pi;
                    ++pj;
                } else if (ci < cj) {
                    ++pi;
                } else {
                    ++pj;
                }
            }

            LLt[(size_t)i * n + j] = s;
        }
    }

    // Compute relative Frobenius error
    double err2 = 0.0;
    double normA2 = 0.0;
    double max_abs_err = 0.0;

    for (size_t k = 0; k < (size_t)n * (size_t)n; ++k) {
        double diff = Adense[k] - LLt[k];
        err2 += diff * diff;
        normA2 += Adense[k] * Adense[k];
        double ad = fabs(diff);
        if (ad > max_abs_err) max_abs_err = ad;
    }

    double rel_frob = (normA2 > 0.0) ? sqrt(err2 / normA2) : 0.0;

    printf("\n");
    printf("====== Numeric Cholesky check ======\n");
    printf("n              = %u\n", n);
    printf("rel |A-LLT| = %.12e\n", rel_frob);
    printf("max |A-LLT|   = %.12e\n", max_abs_err);
    printf("====================================\n");
    printf("\n");
}
#endif

#if !USE_FEM_MATRIX
static void debug_test_solver(const CSRMatrix& A,
                              const CSRMatrix& L,
                              const std::vector<std::vector<ColAdjEntry>>& col_adj)
{
    const uint32_t n = (uint32_t)A.rows;

    std::vector<double> x_true(n);
    std::vector<double> b(n);
    std::vector<double> y(n);
    std::vector<double> x(n);

    // vecteur exact
    for (uint32_t i = 0; i < n; ++i)
        x_true[i] = sin((double)i);

    // b = A x_true
    A.mvp(x_true.data(), b.data());

    // résoudre Ly = b
    forward_solve_lower_csr(L, b.data(), y.data());

    // résoudre LTx = y
    backward_solve_lowerT_with_coladj(L, col_adj, y.data(), x.data());

    // erreur
    double err2 = 0.0;
    double norm2 = 0.0;
    double max_err = 0.0;

    for (uint32_t i = 0; i < n; ++i) {
        double diff = x[i] - x_true[i];
        err2 += diff * diff;
        norm2 += x_true[i] * x_true[i];

        double ad = fabs(diff);
        if (ad > max_err) max_err = ad;
    }

    double rel_err = sqrt(err2 / norm2);

    printf("\n");
    printf("====== Linear solve check ======\n");
    printf("rel |x-x_true| = %.12e\n", rel_err);
    printf("max |x-x_true|   = %.12e\n", max_err);
    printf("================================\n");
    printf("\n");
}
#endif


#if !USE_FEM_MATRIX
void NavierStokesSolver::init_direct_solver(double dt, double nu)
{
    const double alpha = nu * dt;

    if (chol_ready && std::fabs(alpha - chol_alpha) < 1e-12)
        return;

    chol_alpha = alpha;

    // Build A = M + alpha * S 
    build_A_from_M_S(M, S, alpha, A);

    // Symbolic phase
    elimination_tree_from_csr_symmetric(P, etree_parent);
    symbolic_cholesky_row_pattern(P, etree_parent, PL);

    // Numeric phase (Algorithm 5.7 up-looking row-wise)
    cholesky_factorize_rowwise(A, PL, L);

    // Debug numeric phase
    debug_check_cholesky_factorization(A, L);   

    // Build column adjacency for fast L^T solve
    build_col_adjacency_lower(L, col_adj);

    // debug solver
    debug_test_solver(A, L, col_adj);

    chol_ready = true;
}
#endif


void NavierStokesSolver::time_step(double dt, double nu)
{
    size_t iter_stream = compute_stream_function();

    /**********************************************************************
     * Solve the system :
     *
     *  (M + \nu * dt * S)omega(t+dt) = M * omega(t) + dt * T(Omega,Psi)(t)
     *
     *********************************************************************/

#if !USE_FEM_MATRIX

#ifndef NS_USE_DIRECT_SOLVER
#define NS_USE_DIRECT_SOLVER 1
#endif

#if NS_USE_DIRECT_SOLVER
    init_direct_solver(dt, nu);

    double *rhs = p.data;
    double *omega_vec = omega.data;

    compute_transport(rhs);
    blas_axpby(1.0, Momega.data, dt, rhs, N);

    forward_solve_lower_csr(L, rhs, ysolve.data);
    backward_solve_lowerT_with_coladj(L, col_adj, ysolve.data, omega_vec);
#else
    // ---- ancien code CG pour omega ----
    double *residual = r.data;
    double *direction = p.data;
    double *Adir = Ap.data;
    double *omega_vec = omega.data;
    double *Momega_vec = Momega.data;

    compute_transport(direction);
    blas_axpby(1.0, Momega_vec, dt, direction, N);

    double rhs_norm2 = blas_dot(direction, direction, N);

    S.mvp(omega_vec, residual);
    blas_axpby(1.0, Momega_vec, nu * dt, residual, N);
    blas_axpby(1.0, direction, -1.0, residual, N);

    blas_copy(residual, direction, N);

    double r_norm2 = blas_dot(residual, residual, N);
    double rel_error = sqrt(r_norm2 / rhs_norm2);

    size_t iter = 0;

    while (rel_error > tol && iter < iter_max) {
        S.mvp(direction, Adir);
        M.mvp(direction, Momega_vec);
        blas_axpby(1.0, Momega_vec, nu * dt, Adir, N);

        double denom = blas_dot(direction, Adir, N);
        double alpha = r_norm2 / denom;

        blas_axpy(alpha, direction, omega_vec, N);
        blas_axpy(-alpha, Adir, residual, N);

        double new_r_norm2 = blas_dot(residual, residual, N);
        rel_error = sqrt(new_r_norm2 / rhs_norm2);

        if (rel_error <= tol)
            break;

        double beta = new_r_norm2 / r_norm2;
        blas_axpby(1.0, residual, beta, direction, N);

        r_norm2 = new_r_norm2;
        iter++;
    }

    (void)iter;
#endif

#endif

    set_zero_mean(omega.data);
    t += dt;
    (void)iter_stream;

}