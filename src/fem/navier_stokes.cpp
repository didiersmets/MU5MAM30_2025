#include <assert.h>
#include <stddef.h>
#include <stdint.h>

#include "navier_stokes.h"

#include "conjugate_gradient.h"

#include "P1.h"
#include "P2.h"
#include "tiny_blas.h"

NavierStokesSolver::NavierStokesSolver(Mesh &m, const bool &use_fem_P2)
	: m(m), omega(0), Momega(0), psi(0), r(0), p(0), Ap(0), use_fem_P2(use_fem_P2)
{
#if USE_FEM_MATRIX
	if (use_fem_P2)
		throw std::runtime_error("P2 Lagrange elements are not available for FEM matrix framework");
	else
	{
		build_P1_mass_matrix(m, M);
		build_P1_stiffness_matrix(m, S);
	}
#else
	if (!use_fem_P2)
	{
		build_P1_CSRPattern(m, P);
		build_P1_mass_matrix(m, P, M);
		build_P1_stiffness_matrix(m, P, S);
	}
	else
	{
		m.build_edges();
		build_P2_CSRPattern(m, P);
		build_P2_mass_matrix(m, P, M);
		build_P2_stiffness_matrix(m, P, S);
	}
#endif
	rel_error = (double *)malloc(sizeof(double));
	t = 0;
	vol = M.sum();
	inited = false;

	N = M.cols;
	omega.resize(N);
	Momega.resize(N);
	psi.resize(N);
	r.resize(N);
	p.resize(N);
	Ap.resize(N);
}

void NavierStokesSolver::set_zero_mean(double *V)
{
	/* We use the formula : \bar v = (\sum_{i, j} V_i * M_{ij})  / vol
	Then we set V <- V - \bar v */
	double *TEMP = (double *)malloc(N * sizeof(double));
	M.mvp(V, TEMP);
	double sum = blas_sum_in_place(TEMP, N);
	for (size_t i = 0; i < N; i++)
		V[i] -= sum / vol;
}

void NavierStokesSolver::compute_transport(double *T)
{
	/* We use the formula :
	\forall j \in I, T[j] = \sum_{i, k} \Omega_i * \Psi_k \int_{\Omega} \phi_i * (\nabla^T \phi_k . \nabla \phi_j) */
	memset(T, 0, N * sizeof(double));

	double *OMEGA = omega.data;
	double *PSI = psi.data;

	size_t nt = m.triangle_count();
	for (size_t tri = 0; tri < nt; tri++)
	{
		uint32_t a = m.indices[3 * tri];
		uint32_t b = m.indices[3 * tri + 1];
		uint32_t c = m.indices[3 * tri + 2];

		/* Compute T[a] contribution */
		T[a] += (OMEGA[a] * (PSI[c] - PSI[b])) / 6;
		T[a] += (OMEGA[b] * (PSI[c] - PSI[b])) / 6;
		T[a] += (OMEGA[c] * (PSI[c] - PSI[b])) / 6;

		/* Compute T[b] contribution */
		T[b] += (OMEGA[a] * (PSI[a] - PSI[c])) / 6;
		T[b] += (OMEGA[b] * (PSI[a] - PSI[c])) / 6;
		T[b] += (OMEGA[c] * (PSI[a] - PSI[c])) / 6;

		/* Compute T[c] contribution */
		T[c] += (OMEGA[a] * (PSI[b] - PSI[a])) / 6;
		T[c] += (OMEGA[b] * (PSI[b] - PSI[a])) / 6;
		T[c] += (OMEGA[c] * (PSI[b] - PSI[a])) / 6;
	}
}

size_t NavierStokesSolver::compute_stream_function()
{
	size_t iter = 0;

	/**********************************************************************
	 * Solve the system :
	 *
	 *  S * \Psi(t) = M * \Omega(t)
	 *
	 *********************************************************************/

	double *PSI = psi.data;
	double *MOMEGA = Momega.data;
	double *OMEGA = omega.data;
	double *R = r.data;
	double *P = p.data;
	double *AP = Ap.data;

	M.mvp(OMEGA, MOMEGA);

	/* Solve using CG */
	iter = conjugate_gradient_solve(S, MOMEGA, PSI, R, P, AP, rel_error, tol, iter_max, false);

	return iter;
}

void NavierStokesSolver::time_step(double dt, double nu)
{
	double *T = (double *)malloc(N * sizeof(double));
	compute_stream_function();
	compute_transport(T);

	/**********************************************************************
	 * Solve the system :
	 *
	 *  (M + \nu * dt * S) * omega(t+dt) = M * omega(t) + dt * T(Omega,Psi)(t)
	 *
	 *********************************************************************/

	double *MOMEGA = Momega.data;
	double *OMEGA = omega.data;
	double *R = r.data;
	double *P = p.data;
	double *AP = Ap.data;

	double *RHS = (double *)malloc(N * sizeof(double));
	;
	blas_copy(MOMEGA, RHS, N);
	blas_axpy(dt, T, RHS, N);

	/* We directly solve the CG here to avoid creating the matrix M + \nu * dt * S */
	double *TEMP = (double *)malloc(N * sizeof(double));

	b2 = blas_dot(RHS, RHS, N);

	/* r_0 = RHS - M * omega(t) - \nu * dt * S * omega(t) */
	M.mvp(OMEGA, R);
	S.mvp(OMEGA, TEMP);
	blas_axpby(1, RHS, -1., R, N);
	blas_axpby(-nu * dt, TEMP, 1., R, N);

	/* p_0 = r_0 */
	blas_copy(R, P, N);

	r2 = blas_dot(R, R, N);
	*rel_error = sqrt(r2 / b2);

	size_t iter = 0;
	while ((iter < iter_max) && (*rel_error > tol))
	{
		/* Computation of (M + \nu * dt * S) * p_n */
		M.mvp(P, AP);
		S.mvp(P, TEMP);
		blas_axpby(nu * dt, TEMP, 1., AP, N);

		/* Computation of \alpha_n */
		double p_A2 = blas_dot(P, AP, N);
		double alpha = r2 / p_A2;

		/* Computation of omega_{n+1} */
		blas_axpby(alpha, P, 1., OMEGA, N);

		/* Computation of r_{n+1} */
		blas_axpby(-alpha, AP, 1., R, N);

		/* Computation of \beta_{n+1} */
		double new_r2 = blas_dot(R, R, N);
		double beta = new_r2 / r2;

		/* Computation of p_{n+1} */
		blas_axpby(1., R, beta, P, N);

		r2 = new_r2;
		*rel_error = sqrt(r2 / b2);
		iter++;
	}

	set_zero_mean(OMEGA);

	t += dt;
}
