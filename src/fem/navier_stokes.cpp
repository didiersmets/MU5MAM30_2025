#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <vector>

#include "navier_stokes.h"

#include "P1.h"
#include "tiny_blas.h"

NavierStokesSolver::NavierStokesSolver(const Mesh &m)
	: m(m)
	, N(m.vertex_count())
	, omega(N)
	, Momega(N)
	, psi(N)
	, r(N)
	, p(N)
	, Ap(N)
{
#if USE_FEM_MATRIX
	build_P1_mass_matrix(m, M);
	build_P1_stiffness_matrix(m, S);
#else
	build_P1_CSRPattern(m, P);
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

void NavierStokesSolver::time_step(double dt, double nu)
{
	size_t iter_stream = compute_stream_function();

	/**********************************************************************
	 * Solve the system :
	 *
	 *  (M + \nu * dt * S)omega(t+dt) = M * omega(t) + dt * T(Omega,Psi)(t)
	 *
	 *********************************************************************/

	/* Your implementation goes here */

    double *residual = r.data;
    double *direction = p.data;
    double *Adir = Ap.data;
    double *omega_vec = omega.data;
    double *Momega_vec = Momega.data;

    // Construction du second membre

    compute_transport(direction);                 
    blas_axpby(1.0, Momega_vec, dt, direction, N); 

    double rhs_norm2 = blas_dot(direction, direction, N);

    // Résidu initial 

    S.mvp(omega_vec, residual);                   
    blas_axpby(1.0, Momega_vec, nu * dt, residual, N); 
    blas_axpby(1.0, direction, -1.0, residual, N);     

    blas_copy(residual, direction, N);            // p = r

    double r_norm2 = blas_dot(residual, residual, N);
    double rel_error = sqrt(r_norm2 / rhs_norm2);

    size_t iter = 0;

    // Boucle CG pour résoudre A ω^{n+1} = b
    while (rel_error > tol && iter < iter_max) {

        S.mvp(direction, Adir);                 
        M.mvp(direction, Momega_vec);             
        blas_axpby(1.0, Momega_vec, nu * dt, Adir, N); 

        double denom = blas_dot(direction, Adir, N);
        double alpha = r_norm2 / denom;

        // ω = ω + alpha * p
        blas_axpy(alpha, direction, omega_vec, N);

        // r = r - alpha * A p
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

    set_zero_mean(omega.data);

    t += dt;

    (void)iter_stream;
}
